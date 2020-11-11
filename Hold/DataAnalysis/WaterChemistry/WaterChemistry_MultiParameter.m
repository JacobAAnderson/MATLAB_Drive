% Make Water Property Map
clear all;
close all;
clc

global geoImage;
global geoData;
global shape;


%% Check for existing input preferances
savedInputPath = fullfile( userpath,'Script_Inputs');
savedInputs    = fullfile( userpath,'Script_Inputs',[mfilename,'_suggestedFilePaths.mat']);

if exist(savedInputs,'file')
    load(savedInputs);
    
elseif ~exist(savedInputPath,'dir')
    mkdir(savedInputPath);
    suggestions ={pwd,pwd,{'1.0',' -60', '40', '2'},pwd};
    
else
    suggestions ={pwd,pwd,{'1.0',' -60', '40', '2'},pwd};
end


%% Get Data

% Get data to be processed
disp('Select Data File')
try
    [DataFileName,DataFilePath] = uigetfile({'*.mat';'*.log'},'Select MATLAB Data File',suggestions{1});
catch
    [DataFileName,DataFilePath] = uigetfile({'*.mat';'*.log'},'Select MATLAB Data File');
end

if DataFileName == 0                % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                          % End script if there isn't an input for it to use
end


% Get Geo-tiff -----------------------------------------------------
disp('Select Geotif')
try
    [GeotiffName,GeotiffFilePath] = uigetfile({'*.tif';'*.tiff'},'Select Geo-Tiff',suggestions{2});
catch
    [GeotiffName,GeotiffFilePath] = uigetfile({'*.tif';'*.tiff'},'Select Geo-Tiff');
end


if GeotiffName == 0
    GeotiffFilePath = suggestions{2};      % Retain pervious suggested file path if the input box was canceled
    
else
    [geoImage, geoData] = geotiffread([GeotiffFilePath,GeotiffName]);
    
    if class(geoData) ~= 'map.rasterref.GeographicCellsReference'           % Check that geoData is in the right format for 'geoshow'
        
        geoData = georasterref('RasterSize',size(geoImage), ...             % Reformat geoData for 'geoshow'
            'RasterInterpretation', geoData.RasterInterpretation, ...
            'LongitudeLimits',geoData.XWorldLimits, ...
            'LatitudeLimits',geoData.YWorldLimits, ...
            'ColumnsStartFrom',geoData.ColumnsStartFrom,...
            'RowsStartFrom',geoData.RowsStartFrom);
        
    end
end


% get shape for the filter area
disp('Select Shape File')
try
    [shapeName,shapeFilePath] = uigetfile('*.txt','Select Shape File',suggestions{2});
catch
    [shapeName,shapeFilePath] = uigetfile('*.txt','Select Shape File');
end

if shapeName ~= 0
    shape = ImportWayPoints(shapeFilePath, shapeName,'\n');
else
    shape = 0;
end


% Get mapping specifications
disp('Set EKF specifications')
fields = {'Grid cell width','Floor of the map','Window size of EKF','Sigmas'};
try
    answers = inputdlg(fields, 'Input Specs...',[1 20; 1 20; 1 20; 1 20], suggestions{3} );
catch
    answers = inputdlg(fields, 'Input Specs...',[1 20; 1 20; 1 20; 1 20], {'1.0',' -60', '40', '2'} );
end

if length(answers) < 1
    answers = suggestions{3};
end


% Save the Inputs for suggestions next time the script is called
% This happens twice, once here incase the script is treminated before compleation and again at the end of the scrip after the final dialog box
suggestions = {DataFilePath,GeotiffFilePath,answers,suggestions{4}};
save(savedInputs,'suggestions');


%% Load and Sort the Data

dataFile = fullfile(DataFilePath,DataFileName);

if contains(DataFileName,'ecomapper','IgnoreCase',true)
    
    disp('Loading data')
    disp('This may take a few minuts')
    load( dataFile );
    
    % Extract data from data structure
    vehicle = cat(1,data.vehicle);                    % Latitude, Longitude, Depth from surface
    Data    = cat(1,data.wqData);                     % Water chemistry
    Date    = cat(1,data.date);                       % Date and time
    header  = strsplit(data(4).Header,{':',','});     % Header for the water chemistry matrix
    header  = header(2:end);                          % Eliminate the header label from the header
    
    clear('data')                                     % Free up some space after the data has been extracted from the storage stucture
    
    % Sort data chronologicaly
    disp('Sorting the data chronologically')
    [Date, index] = sortrows(Date,'ascend');
    
    vehicle = vehicle(index,:);
    Data    = Data(index,:);
    
    % Get stats. on the data being processed
    dataMax  = max(Data,[],1);
    dataMean = mean(Data,1);
    dataSTD  = std(Data,1);
    dataMin  = min(Data,[],1);
    
    % Display a summery of the table
    f = figure('Name','Data Summery','NumberTitle','off');
    t = uitable(f,'Data',[dataMax; dataMean; dataMin; dataSTD], 'Units','Normalized','Position',[0 0 1 1]);
    t.RowName = {'Data Max', 'Data Mean', 'Data Min', 'Data STD'};
    t.ColumnName = header;
    
    % Save the sorded data for reuses
    disp('Save the Sorted Data')
    
    stats = [dataMax, dataMean, dataSTD, dataMin ];
    uisave({'vehicle', 'Data','Date','header','stats','f'},[num2str(yyyymmdd(datetime('today'))),'_ChronologicalySortedData']);
    
    
elseif contains(DataFileName,'ChronologicalySorted','IgnoreCase',true)
    disp('Loading Previously Sorted data')
    load( dataFile );
    
    % Get stats. on the data being processed
    dataMax  = stats(1);
    dataMean = stats(2);
    dataSTD  = stats(3);
    dataMin  = stats(4);
end


%% Choose how to save outputs
disp('Where to Save Animation')
try
    [videoName,videoPath] = uiputfile('*.mp4','Save Aninations As',[suggestions{4},num2str(yyyymmdd(datetime('today'))),'_Movie.mp4']);
catch
    [videoName,videoPath] = uiputfile('*.mp4','Save Animations As',[pwd,num2str(yyyymmdd(datetime('today'))),'_',name,'_Movie.mp4']);
end

if videoName == 0
    disp('Animations Have Not Been Saved')
    return
    
else
    % Save the Inputs for suggestions next time the script is called
    suggestions = {DataFilePath,GeotiffFilePath,answers,videoPath};
    save(savedInputs,'suggestions');
    
end


%% Cycle through data in one hour Increments
for parameter = 1: length(header)
    
    depth = vehicle(:,3)<= 0.25;
    
    % Pull out the values that we are interesedted in
    
    % % All depth
    % latitude  = vehicle(:,1);
    % longitude = vehicle(:,2);
    % Data      = Data(:,parameter);
    
    % % Filter by depth
    latitude  = vehicle(depth,1);
    longitude = vehicle(depth,2);
    Data      = Data(depth,parameter);
    
    units     = header(parameter);
    
    ii = 1;
    count = 1;
    startNum = Date(1);
    
    h = waitbar(0,'Processing Data');
    for index = 1:length(Data)
        
        if abs(Date(index) - startNum) >= hours(2) || index == length(Data)
            
            fig = MakeMape(hr_lon', hr_lat', hr_data', answers, units, startNum, dataMax(parameter), dataMin(parameter),dataMean(parameter), dataSTD(parameter));
            Frames(ii) = getframe(fig);
            close(fig);
            
            startNum = Date(index);
            count = 1;
            ii = ii+1;
            
            waitbar(index / length(Data));
            
        end
        
        hr_lat(count)  = latitude(index);
        hr_lon(count)  = longitude(index);
        hr_data(count) = Data(index);
        count = count+1;
        
    end
    
    try
        close(h);
    end
    
    
    % Play sound to indicate that thge script has finnished
    sound_ = load( 'gong');
    sound(sound_.y, sound_.Fs)
    
    
    % Alow user to view the plot before being prompted to save the data
    disp('MATLAB is paused!!!')
    disp('Press any button to continue')
    pause
    
    
    % Play animation of the figures
    movFig = figure('Name','Water Chemistry Animation','NumberTitle','off');
    movAxes = axes(movFig,'Units', 'normalized','Position',[0.04 0 0.7 0.7] );
    axis off
    movie(movAxes,Frames,1,1);
    
    
    % Save animation as a movie
    try
        name = strsplit(units{:},{'['});   % Get just the name of the data peramiter
    catch
        name = strsplit(units{:},{'+'});   % Get just the name of the data peramiter
    end
    
    name = name{1};                 % Only keep the firat cell
    name = name(~isspace(name));    % Eliminate any white spaces in the string
    
    % Get file path to save the video
    
    
    videoPathName = fullfile(videoPath,videoName);  % make full file path for the video
    
    % Create a video from the animation
    v = VideoWriter(videoPathName,'MPEG-4');
    v.FrameRate = 1;
    open(v);
    writeVideo(v,Frames);
    close(v)
    
    % Save the frames for the move for individual access later
    [framesFilepath,framesName,framesExt] = fileparts(videoPathName);
    saveFrames = fullfile(framesFilepath,[framesName,'_Frames.mat']);
    save(saveFrames,'Frames')
    
end






%% Functions ======================================================================================================================================================
% =================================================================================================================================================================
% =================================================================================================================================================================


%% Make Water Chemistry Map
function fig = MakeMape(longitude, latitude, Data, answers, units, startNum, dataMax, dataMin, dataMean, dataSTD )

global geoImage;
global geoData;
global shape;


% disp('Making water chemistry Map')
res       = str2double(answers{1}); % Grid cell width
baseFloor = str2double(answers{2}); % Floor of the map --> This will defalt to the minimum data value if the argument is not passed into the Bathymetry function.
window    = str2double(answers{3}); % Window size of EKF
sig       = str2double(answers{4}); % Sigma threshold for low pass filter

[Bathymetry, LON,  LAT ] = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', window, 'sigma', sig, 'shape',shape); % ,'Floor', baseFloor );


% Get the variable name out of Units
try
    name = strsplit(units{:},{'['});   % Get just the name of the data peramiter
catch
    name = strsplit(units{:},{'+'});   % Get just the name of the data peramiter
end

name = name{1};                 % Only keep the firat cell
name = name(~isspace(name));    % Eliminate any white spaces in the string

figName = [char(datetime(startNum,'ConvertFrom','datenum','Format','MMMM dd, yyyy HH:mm:ss')),' ',name];



% Plot Water chemistry ----------------------------------------------------
try
    fig = figure('Name',figName,'NumberTitle','off','Visible', 'off');
    geoPlot = geoshow(geoImage,geoData);
    
    hold on
    
catch
    fig = figure('Name',figName,'NumberTitle','off','Visible', 'off');
end



%-- Filled Contours -------------------------------------------------------
[c,h] = contourf(LON, LAT,Bathymetry, 'k');
clabel(c,h,'FontSize',8)

if abs(dataSTD) > 10
    caxis([dataMin, round(dataMean+3*dataSTD)]);
else
    caxis([(dataMin-2), (dataMax+2)]);
end

cbar = colorbar;
hold off
title(figName)
xlabel('Easting [deg]')
ylabel('Norting [deg]')
zlabel('Depth [meters]')
title(cbar,units)



% Save figures
% figName = [char(datetime(startNum,'ConvertFrom','datenum','Format','yyyyMMdd_HHmmss')),'_',name];
% savePath = fullfile('C:\Users\jaanderson.FORTLEWIS\_OutPutts',figName);
% savefig(savePath);

end



