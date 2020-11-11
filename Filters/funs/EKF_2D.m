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
function [meanGrid, varargout] = EKF_2D(dataX, dataY, dataRaw, varargin)

% Initalize parameters that can be change by the varargin --------------------------------------------------------------
gridCellResolution = 0.1;
window = 20;
baseFloor = min(dataRaw(:));
std_theshold = 2;
sigma_slope = 0.1204;
sigma_offset = 0.6142;
shape = 0;

% Get Additional Input peramiters ---------------------------------------------------------------------------------------
if  nargin > 1
    
    for i=1:length(varargin)
        try
            
            switch lower(cell2mat(varargin(i)))
                
                case 'window'                               % Get Window size
                    window = varargin{i+1};
                    
                case 'gridcell'                             % Get grid cell resolution
                    gridCellResolution = varargin{i+1};
                    
                case'floor'                                 % Get baes floor level
                    baseFloor = varargin{i+1};
                    
                case 'sigma'                                % Get sigma threshhold for filtering outliers
                    std_theshold = varargin{i+1};
                    
                case 'slope'                                % Get the slope for the covariance calculataion
                    sigma_slope = varargin{i+1};
                    
                case 'offset'                               % Get the offset for the covariance calculataion
                    sigma_offset = varargin{i+1};
                    
                case 'shape'                                % Get outline for the filter area
                    shape = varargin{i+1};
                    
                otherwise
                    continue
            end
            
        catch
            continue
            
        end
    end
end


% Filter out outliers from altitude measurements -----------------------------------------------------------------------------------

dataFiltered(1) = dataRaw(1);                             % Create a second Data array for filtered data values
dataFiltered(length(dataRaw)) = dataRaw(length(dataRaw));

r = 5;                                          % Number of points, on either side of the data point being evaluated, used to determin the STD and mean of the data point

for j = (r+1): length(dataRaw)-(r+1)               % Cycle through the data set.
    % Reassing values that are "std_theshold" greater than the STD to be the average value
    s = std(dataRaw(j-r:j+r));
    u = mean(dataRaw(j-r:j+r));
    
    if abs(dataRaw(j)-u) > std_theshold*s
        dataFiltered(j)= u;
    else
        dataFiltered(j) = dataRaw(j);
    end
end



% Create boundary map ------------------------------------------------------------------------------------

% Finde the origen of the data set
origin(1) = min(dataX);
origin(2) = min(dataY);

% Deteruming the oposite corner of the data grid
Lmax(1) = max(dataX);
Lmax(2) = max(dataY);


% Get the shape of the area coverd by the data as a ploygon
if shape == 0                       % Boundary polygon was not supplied by the user
    
    in = boundary(dataX,dataY);                                             % Find the outer edge of the data points
    edgePoints = [dataX(in),dataY(in); 0,0];                                % Close the polygon
    edgePoints(end,1) = edgePoints(1,1);
    edgePoints(end,2) = edgePoints(1,2);
    
    edgePoints = expandPolygon([edgePoints(:,1), edgePoints(:,2)], window* gridCellResolution/2 );  % Expand th epolygon by half the windo size
    edgePoints = edgePoints{:};
    
    if any(isinf(edgePoints(:))) || any(isnan(edgePoints(:)))               % Filter out any unwanted numbers --> inf, nan
        
        edgePoints(any(isinf(edgePoints),1),:) = [];
        edgePoints(any(isinf(edgePoints),2),:) = [];
        edgePoints(any(isnan(edgePoints),1),:) = [];
        edgePoints(any(isnan(edgePoints),2),:) = [];
        
    end
    
    if edgePoints(1,1) ~= edgePoints(end,1) || edgePoints(1,2) ~= edgePoints(end,2)     % Close the polygon
        edgePoints = [ edgePoints(:,1), edgePoints(:,2)
            edgePoints(1,1), edgePoints(1,2)];
        
    end
    
%     % Plot the data boundary to check the polygon expansion
%     figure('Name','EKF_2D Expaneded Boundary')
%     hold on
%     plot(edgePoints(:,1), edgePoints(:,2),'r')
%     plot(dataX(in),dataY(in),'b')
%     legend('Expaneded Boundary','Data Boundary')
%     hold off
    
    
else                                 % Boundary polygon was supplied by the user
    edgePoints = shape;
    if edgePoints(1,1) > 0
        edgePoints = fliplr(edgePoints);
    end
end


% Finde the origen of the interpolation area
origin(1) = min(edgePoints(:,1));
origin(2) = min(edgePoints(:,2));

% Deteruming the oposite corner of the interpolation area
Lmax(1) = max(edgePoints(:,1));
Lmax(2) = max(edgePoints(:,2));

% Size of the map
dims = Lmax - origin;
n = ceil(dims(1)/gridCellResolution);   % Coulmns
m = ceil(dims(2)/gridCellResolution);   % Rows

% Create Area mesh
x_ax = linspace(origin(1),Lmax(1), n);
y_ax = linspace(origin(2),Lmax(2), m);
[xGrid, yGrid] = meshgrid(x_ax, y_ax);

% Create a map to indicate which elements of the mean and variance matices should be updated
dataArea    = inpolygon(xGrid,yGrid,edgePoints(:,1),edgePoints(:,2));
updatedArea = false(size(dataArea));


% Initialize Bathymetry map -----------------------------------------------------------------------------------------------
meanGrid = baseFloor * ones(m,n);
variance = zeros(m,n);


% EKF ------------------------------------------------------------------------------------------------------------------------
% loop over all elements of auv path L, to fuse the altimeter measurements with each cell within win cells of auv position
for t = 1:length(dataFiltered)
    
    indicies = getIndices([dataX(t),dataY(t)], origin, Lmax, [n, m]);              % Get index of lat,long of t^th element of L
    j = indicies(1);
    i = indicies(2);
    
    % update all cells with the window of cell i,j
    for row = max(1,i-window): min(m,i+window)
        
        for col = max(1,j-window): min(n,j+window)
            
            % Check if we care to update this part of the map, i.e. it
            % is within boundaries of our edgepoints
            if dataArea(row,col) == 1
                
                updatedArea(row,col) = true;
                % Use a 1D Kalman Filter for each cell to estimate depth to seafloor.
                % Variance of normally distributed random varialbe representing each measurement is a
                % function of planar distance between AUV position and cell position.
                dist2 = 0.1+gridCellResolution*((row-i)^2+(col-j)^2)^0.5;
                
                % Check if this is the first measurement associated
                % with this cell
                if meanGrid(row,col) == baseFloor
                    meanGrid(row,col) = dataFiltered(t);
                    variance(row,col) = (dist2*sigma_slope+sigma_offset)^2;
                    
                    % Fuse new measurements for previously
                else
                    var=(dist2*sigma_slope+sigma_offset)^2;
                    K = variance(row,col)/(variance(row,col) + var);
                    meanGrid(row,col) = meanGrid(row,col)+ K*( dataFiltered(t) - meanGrid(row,col)) ;
                    variance(row,col) = variance(row,col)- K*variance(row,col);
                end
            end
        end
    end
    
end

meanGrid(~updatedArea) = nan;    % Clear any readings outside of the data area



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
function indicies = getIndices(X,  origin, uperRight, dims)

diffs = (X - origin) ./ (uperRight - origin);

indicies = floor( dims .* diffs )+1;

end

