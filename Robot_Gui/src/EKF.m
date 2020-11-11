%% EKF_2D is a one dimentional Extended Kalman filter applied tin two dimentions.

% The filter will first loop through the data and remove outliers that are more than n sigmas away from the aveage value.
% n is can be specified as an input into the function. Its default values is 2

% The filter will then apply an Extended kalman filter to a two dimentioan, rectangular, window.
% The size of the window can be specified as an input but the default size is 20 units on either side of the data point being evaluated

% The filter is curently set up for taking longitudes and latitudes as x y values

% INPUTS:
%   x    --> longitudinal values
%   y    --> latitudinal values
%   Data --> Data to be used of interpolation

% Optional Inputs:
%   gridCellResolution  --> Size of the data grid: value of 1 means each cell will be one unit, value of 2 means each cell will be 2 units
%   window              --> Size of the window used for data intermpolation
%   baseFloor           --> Number used to initialized meanGrid, defaul is the minimunvalu of the data set
%   std_theshold        --> Vumber of sigmas to threshold on when eliminating outliers
%   sigma_slope         -->
%   sigma_offset        -->

% OUTPUTS
%   meanGrid --> interpoatoed data values

% Optional Outputs:
%   xGrid    --> Mesh of the x coordiantes associated with the interpolated dataGrid
%   yGrid    --> Mesh of the y coordiantes associated with the interpolated dataGrid
%   variance --> Variance of the interploalated data points

%__Syntax_______________________________________________________________________________________________________________________________________
%
%                            meanGrid = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', window, 'sigma', sig, 'Floor', baseFloor );
%               [meanGrid, variance ] = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', window, 'sigma', sig, 'Floor', baseFloor );
%          [meanGrid, xGrid,  yGrid ] = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', window, 'sigma', sig, 'Floor', baseFloor );
% [meanGrid, xGrid,  yGrid, variance] = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', window, 'sigma', sig, 'Floor', baseFloor );
%_______________________________________________________________________________________________________________________________________________



%%
function [meanGrid, varargout] = EKF(x, y, Data, varargin)
global ekfProperties;
% Initalize parameters that can be change by the varargin --------------------------------------------------------------

window = 20;
baseFloor = min(Data(:));
std_theshold = 2;
sigma_slope = 0.1204;
sigma_offset = 0.6142;


% Get Additional Input peramiters ---------------------------------------------------------------------------------------
if  nargin > 1
    
    for i=1:length(varargin)
        try
            
            switch lower(cell2mat(varargin(i)))
                
                case 'window'                               % Get Window size
                    window = varargin{i+1};
                    
                   
                case'floor'                                 % Get baes floor level
                    baseFloor = varargin{i+1};
                    
                case 'sigma'                                % Get sigma threshhold for filtering outliers
                    std_theshold = varargin{i+1};
                    
                case 'slope'                                % Get the slope for the covariance calculataion
                    sigma_slope = varargin{i+1};
                    
                case 'offset'                               % Get the offset for the covariance calculataion
                    sigma_offset = varargin{i+1};
                    
                   
                otherwise
                    continue 
            end
            
        catch
            continue

        end
    end
end


% Filter out outliers from altitude measurements -----------------------------------------------------------------------------------

Data2(1) = Data(1);                             % Create a second Data array for filtered data values
Data2(length(Data)) = Data(length(Data));

pOut = 0;                                       % Keep track of the number of points that are filtered out
r = 5;                                          % Number of points, on either side of the data point being evaluated, used to determin the STD and mean of the data point

for j = (r+1): length(Data)-(r+1)               % Cycle through the data set.
    % Reassing values that are "std_theshold" greater than the STD to be the average value
    s = std(Data(j-r:j+r));
    u = mean(Data(j-r:j+r));
    
    if abs(Data(j)-u) > std_theshold*s
        pOut = pOut+1;
        Data2(j)= u;
    else
        Data2(j) = Data(j);
    end
end


% Initialize Bathymetry map -----------------------------------------------------------------------------------------------
Origin             = ekfProperties{1};
gridCellResolution = ekfProperties{2};
dataArea           = ekfProperties{3};

[m,n] = size(dataArea);
meanGrid = baseFloor*ones(m,n);
variance = zeros(m,n);
updatedArea = false(m,n);

% EKF ------------------------------------------------------------------------------------------------------------------------
% loop over all elements of auv path L, to fuse the altimeter measurements with each cell within win cells of auv position
% h = waitbar(0,'Processing Data','Name','EKF_2D');
for t = 1:length(Data2)
    
    [x_pos,y_pos] = getXYCoordinates([x(t),y(t)], Origin);                       % Get x,y position of the data point in meters
    [i,j] = getIndices([x_pos y_pos], gridCellResolution);              % Get index of lat,long of t^th element of L
    
    
    % update all cells with the window of cell i,j
    for k = max(1,i-window):1:min(m,i+window)
        
        for l = max(1,j-window):1:min(n,j+window)
            
            % Check if we care to update this part of the map, i.e. it
            % is within boundaries of our edgepoints
            if dataArea(k,l) == 1
                
                updatedArea(k,l) = true;
                % Use a 1D Kalman Filter for each cell to estimate depth to seafloor. 
                % Variance of normally distributed random varialbe representing each measurement is a
                % function of planar distance between AUV position and cell position.
                dist2 = 0.1+gridCellResolution*((k-i)^2+(l-j)^2)^0.5;
                
                % Check if this is the first measurement associated
                % with this cell
                if meanGrid(k,l) == baseFloor
                    meanGrid(k,l) = Data2(t);
                    variance(k,l) = (dist2*sigma_slope+sigma_offset)^2;
                    
                    % Fuse new measurements for previously
                else
                    var=(dist2*sigma_slope+sigma_offset)^2;
                    K = variance(k,l)/(variance(k,l) + var);
                    meanGrid(k,l) = meanGrid(k,l)+ K*( Data2(t) - meanGrid(k,l)) ;
                    variance(k,l) = variance(k,l)- K*variance(k,l);
                end
            end
        end
    end
    
    % waitbar( t / length(Data2))
    
end

meanGrid(~updatedArea) = nan;    % Clear any readings outside of the data area

% % Close the wait bar
% try
%     close(h)
% catch exception
%     disp(exception.message)
% end

% Assign any additional outputs
if nargout > 1
    switch nargout
        
        case 2
            varargout{1} = variance;
            
        case 3
            varargout{1} = xGrid;
            varargout{2} = yGrid;
            
        case 4
            varargout{1} = xGrid;
            varargout{2} = yGrid;
            varargout{3} = variance;
    end
end

end



%% Determine the row,column indicies of a cell given an x,y location
% X = [ x , y ]
function [row,column] = getIndices(X, gridCellResolution)

row    = floor(X(:,2)./gridCellResolution)+1;
column = floor(X(:,1)./gridCellResolution)+1;
end


%% Converts latitude to an x value, and longitude to a y value
%       L is the 2x1 vector with latitude, longitude
%       L0 is the 2x1 vector with the latitude, longitude of the coordinate frame origin
function [x,y] = getXYCoordinates(L, Origin)

lat  = L(:,2);
long = L(:,1);

[m,~] = size(L);

lat0(1:m,1) = Origin(2);
long0(1:m,1) = Origin(1);

x = vdist(lat0,long,lat0,long0); % from longitude to x in meters
y = vdist(lat,long0,lat0,long0); % from latitude  to y in meters

x(isnan(x)) = 0;
y(isnan(y)) = 0;

end

