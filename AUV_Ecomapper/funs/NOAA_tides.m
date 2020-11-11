% NOAA Tides
% Jacob Anderson
% 8/15/2019

% Get Tide Data from noaa at a specified staion and time duration

% station     --> NOAA station ID [ doulbe ], find at  https://tidesandcurrents.noaa.gov/products.html
% start       --> Beging of time duration [ datetime ]
% stop        --> End of time duration [ datetime ]
% plotResults --> Optianal input to indicat that the interpolated results should be plotted [ boolean ]



function tide = NOAA_tides( station, start, stop, plotResults)


% Build URLs ---------------------------------------------------------------------------
noaa     = 'https://tidesandcurrents.noaa.gov';
apicall  = '/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL';
duration = sprintf('&begin_date=%d&end_date=%d', yyyymmdd(start), yyyymmdd(stop));
location = sprintf('&station=%d', station);
format   = '&datum=MLLW&time_zone=lst_ldt&units=metric&interval=hilo&format=csv';

url_metaData = sprintf('%s/mdapi/v1.0/webapi/stations/%d.json',noaa, station);
url_data = strcat(noaa, apicall, duration, location, format);


fprintf('\n\n_____Geting_Tide_Data_____________________________________________________\n\n')
fprintf('Requesting Data From : %s\n\n', noaa)

% Make Web Request ----------------------------------------------------------------------------
try
    meta_data = webread( url_metaData );                                    % Request meta-data
    
    location = sprintf('%s. %s', meta_data.stations.name, meta_data.stations.state);
    
    fprintf("Accessing Tide Data For:  %s  at %f, %f\n\n", location, meta_data.stations.lat, meta_data.stations.lng)
    
    data = webread( url_data );                                             % Request data
    
    disp("Tide Data Recived")
    disp(data)
    fprintf("\n\n")
    
catch Error
    
    disp(Error)
    fprintf('\n\n !!No Tide Data Recived!! \n\n')
    tide = [];
    return
end


% Interpolate water levels between hight and low tides --------------------------------------
dataTimes  = data.DateTime;
prediction = data.Prediction;
type       = data.Type;

dT = seconds( dataTimes(end) - dataTimes(1));

time = linspace(dataTimes(1), dataTimes(end), dT);
height = zeros(size(time));

start = 1;
for ii = 1: size(data,1) - 1
    
    A = abs(prediction(ii+1) - prediction(ii))/2;
    dT = seconds( dataTimes(ii+1) - dataTimes(ii));
    
    if type{ii} == 'L'
        theta0 = pi;
        thetaN = 0;
        dY = prediction(ii) + A;
        
    elseif type{ii} == 'H'
        theta0 = 0;
        thetaN = pi;
        dY = prediction(ii) - A;
        
    else
        disp("The Tides are Agains Us Master!!")
    end
    
    theta = linspace(theta0, thetaN, dT);
    
    height(start: start + dT - 1) = A * cos(theta) + dY;
    
    start = start + dT;
    
end

height( numel(time)+1: end) = [];               % Get rid of extra entries


% Put Data in structure -------------------------------------------------------------
tide.DateTime = time;
tide.Height = height;
tide.Location.Name = location;
tide.Location.Lat = meta_data.stations.lat;
tide.Location.Lon = meta_data.stations.lng;
tide.Units = 'Meters';



% Display data ----------------------------------------------------------------------
if nargin > 3 && plotResults
    figure('Name','Tide Information','NumberTitle', 'off')
    plot(time, height)
    xlabel("Date-Time")
    ylabel("Tide Height [m]")
    title(location)
end

fprintf('\n\n__________________________________________________________________________\n\n')

end
