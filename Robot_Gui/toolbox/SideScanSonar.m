% Side Scan Soanr Tool
% Started: 07/2019
% Last update : 08/2019


% ___ Side Scan Sonar Data Structure _____________________________________________________________________________
% utm       = cat(1,Data.position);             % Vehicle location in UTMs [ Easting, Northing, Altitude] 
% utmzone   = cat(1,Data.utmZone);              % UTM-zone: Longitude and Latitude bands 
% heading   = cat(1,Data.heading);              % Vehicle heading [deg]
% speed     = cat(1,Data.speed);                % Vehicle speed   [m/s]
% depth     = cat(1,Data.depth);                % Water depth     [m]
% portSonar = cat(1,Data.portSonar);            % Port Sonar data structure
% starSonar = cat(1,Data.starSonar);            % Starboard Sonar data structure
% portEcho  = cat(2,portSonar.EchoStrength);    % Concatinate the image date from the sonar data structures 
% starEcho  = cat(2,starSonar.EchoStrength);
% _________________________________________________________________________________________________________________


%% Side Scan Sonar GUI
function SideScanSonar(fromPanel)
global SSSWindow                                                           % Side Scan Sonar GUI Window

disp('Side Scan Sonar')

% Make SSS GUI Window ----------------------------------------------------------------------------------------------------------------------------------------
windowName = 'Side Scan Sonar';                                             % GUI Windo Name
panelName  = 'Transformations';                                             % Button Panel name

SSSWindow = GUIwindow(windowName);                                          % Instanciate the GUI Windo
SSSWindow.GUIpanels(panelName,'make','position',[0.8 0 0.2 1],'backgroundcolor',[0.4706  0.6706  0.1882]);  % Create putton panel for primary GUI functions
SSSWindow.Header.BackgroundColor = [0.4706  0.6706  0.6];

% Make buttons in the button panel: 
%{ Buttin type,   Button Name,   Callback function and their inputs,   Callback function for slider listeners }
SSSWindow.GUIbuttons( panelName,'add',{
     'pushbutton', 'Load New Data',   {@LoadSSSData, fromPanel, panelName, 0}, [], 'Load A New Soanr File';                                     %   1
     'pushbutton', 'Batch Data Dir.', {@LoadSSSData, fromPanel, panelName, 1}, [], 'Batch Process Directory of Soanr Files';                    %   2
     'pushbutton', 'Load Features',   @LoadFeatures,                           [], 'Load Previousl Anotated Features';                          %   3
     'space',      ' ',               [],                                      [], '';                                                          % --4
     'pushbutton', 'Mark Up',         {@ImageMask,fromPanel},                  [], 'Create Binary Mask of Interesting Locations on the Image';  %   5
     'pushbutton', 'Move Feature',    @AdjustMasks,                            [], 'Create Binary Mask of Interesting Locations on the Image';  %   6
     'pushbutton', 'Merge Feature',   @MergeFeatures,                          [], 'Create Binary Mask of Interesting Locations on the Image';  %   6
     'space',      ' ',               [],                                      [], '';                                                          % --7
     'pushbutton', 'Remove Feature',  @RemoveFeature,                          [], 'Select Feature to Remove';                                  %   8
     'space',      ' ',               [],                                      [], '';                                                          % --9
     'pushbutton', 'Reset Features',  @ResetFeatures,                          [], 'Delete All Features';                                       %  10
     'space',      ' ',               [],                                      [], '';                                                          % -11
%      'pushbutton', 'Crop Image',      @CropImage,                              [], 'Addjusts The Boundary of the Image and the data';           %  10
%      'pushbutton', 'Flip Image',      @FlipImage,                              [], 'Rotate the Image and data 180';                             %  11
     'pushbutton', 'GetPatch',        @GetRectangularPatch,                    [], 'ROI rectangle to get an image patch';                       %  12
     'space',      ' ',               [],                                      [], '';                                                          % -13
     'pushbutton', 'Export Images',   @ExportImage,                            [], 'Save subimages and binary mask';                            %  14
     'pushbutton', 'Save Features',   @SaveFeatures,                           [], 'Save Features ';                                            %  15
     'pushbutton', 'Save SSS data',   @SaveSSSdata,                            [], 'Save Features ';                                            %  15
     });

% Panel = SSSWindow.Panels(panelName);                                        % Extract the button pane from its mapped container to access the structure
% Panel.button(13).Push.Visible = 'off';   

% Panel = SSSWindow.Panels('Transformations');                                % Extract the button pane from its mapped container to access the structure
% Panel.button(13).Push.Visible = 'on';   

end


%% Load Data
% --- Load SSS Data ---------------
function LoadSSSData(~,~, fromPanel, panelName, dir)
global ui SSSWindow DATA 

SSSWindow.Header.String = 'Select Data Input';                              % Prompt to select data
SSSWindow.mouseCounter = 1;

if dir
    filePath = ui.getDir('SonarData','logdoc');                             % Get folder GUI   
    if isempty(filePath), return, end                                       % End script if there isn't an input for it to use 
else
    filePath = ui.getFile('SonarData','logdoc');                            % GUI to get the name and file path of a file
    if isempty(filePath), return, end                                       % End script if there isn't an input for it to use
end


SSSWindow.Header.String = 'Data Is Loading';                                % Prompt to wait for data to finish loading
set(SSSWindow.Figure, 'pointer', 'watch')
drawnow;

if exist(sprintf('%s_SSS_Data.mat', filePath),'file')
    disp("Loading Previously Saved Instance")
    sss = load(sprintf('%s_SSS_Data.mat', filePath));
    DATA.SSS_data = sss.data;
    
else    
    DATA.SSS_data = SSS_Data(filePath);                                     % Read Logdoc File
end

set(SSSWindow.Figure, 'pointer', 'arrow')
drawnow;

if isempty(DATA.SSS_data)                                                   % Check that data is present
    SSSWindow.Header.String = 'Empty Data File';                            % Prompt to select data
    return                                                                  % End script if there isn't any data to process
end

if dir                                                                      % Add a mouse scroll function to togle throught the SSS files
    SSSWindow.AddMouseWheel( @MouseCallBack, {fromPanel, panelName}, 1, size(DATA.SSS_data.Images,2) );
    
    if ~exist(sprintf('%s_SSS_Data.mat', filePath),'file')
        data = DATA.SSS_data;
        disp("Saving Data As a Matlab Variable")
        save(sprintf('%s_SSS_Data.mat', filePath), 'data', '-v7.3')
    end
end

DATA.EM_data = DATA.EM_data.FilterAltimiter(2);

DATA.SSS_data = DATA.SSS_data.Index2EM_Data(DATA.EM_data.RawData.timeStamp);
DATA.SSS_data = DATA.SSS_data.Add_EM_Alt(DATA.EM_data.Altimeter);

DrawSSSImage(fromPanel)                                                     % Display SSS Image

end


% --- Load Features ---------------
function LoadFeatures(~,~)
global DATA ui

filePath = ui.getFile('SSS_Features','mat');                                % GUI to get the name and file path of a file
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use

DATA.SSS_data = DATA.SSS_data.LoadFeatures(filePath);

end


%% Save Data
% --- Save SSSdata object
function SaveSSSdata(~,~)
global DATA ui

filePath = ui.saveFile('SSS_Data','mat');                                % GUI to get the name and file path of a file
if isempty(filePath), return, end                                             % End script if there isn't an input for it to use

data = DATA.SSS_data;
disp("Saving Data As a Matlab Variable")
save(filePath, 'data', '-v7.3')

end


% ---Save Features -------------
function SaveFeatures(~,~)
global DATA ui

filePath = ui.saveFile('SonarFeatures','mat');                                % GUI to get the name and file path of a file
if isempty(filePath), return, end                                             % End script if there isn't an input for it to use

DATA.SSS_data.SaveFeatures(filePath);
end


% --- Export Subimages ----------
function ExportImage(~,~)

global DATA ui

% Save Data
disp("Choose Loaction to Save Images and Binary Masks")

filePath = ui.getDir('SideScan_FeatureImages','jpeg');                % GUI to get file path to save file
if isempty(filePath), return, end

list = DATA.SSS_data.GetFeatureSheetList;

indx = UI_ChecBox("Sonar Data Sets", list);

disp('Saving Images From:')
disp([list(indx)]')

DATA.SSS_data.ExportImagePatches(filePath, list(indx));

end


%% Add, Delete, maipulate Features
% ---Create Feature ---------------------------
function ImageMask(~,~,fromPanel)
global SSSWindow DATA  
pos = SSSWindow.InteractiveLayer;                                          % Create interactive polygon in the GUI and wait for it to be closes
DATA.SSS_data = DATA.SSS_data.AddFeatue(pos, SSSWindow.mouseCounter);
DrawSSSImage(fromPanel)
end


% --- Addust Feature --------------------------------------------
function AdjustMasks(~,~)
global SSSWindow DATA  

[x, y] = getpts(SSSWindow.Axes);

[m,n] =size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask);

ind = sub2ind([m,n], round(y), round(x));

region = DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask(ind);

BW = false(m,n);

BW(DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask == region) = true;

Boundaries = bwboundaries(BW,'holes'); 
Boundaries = Boundaries{1};

[x,y] = reducem( Boundaries(:,2), Boundaries(:,1),20);

% Find the boundary that corresponds to the Occupancy grid
pos = SSSWindow.InteractiveLayer(x ,y);

DATA.SSS_data = DATA.SSS_data.AdjustFeatue(pos, region, SSSWindow.mouseCounter);

DrawSSSImage

end


% --- Reset Features ---------------------------------
function ResetFeatures(~,~)
global DATA
DATA.SSS_data = DATA.SSS_data.ResetData;
DrawSSSImage
end


% --- Remove Feature ------------------------------------
function RemoveFeature(~,~)
global SSSWindow DATA 

[x, y] = getpts(SSSWindow.Axes);

[m,n] =size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask);

ind = sub2ind([m,n], round(y), round(x));

region = DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask(ind);

DATA.SSS_data = DATA.SSS_data.RemoveFeature(region, SSSWindow.mouseCounter);

DrawSSSImage
end


% --- Merge Features ------------------------------------
function MergeFeatures(~,~)

global SSSWindow DATA  

[x, y] = getpts(SSSWindow.Axes);

[m,n] =size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask);

ind = sub2ind([m,n], round(y), round(x));

region = DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask(ind);

% BW = false(m,n);
% 
% BW(DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask == region) = true;
% 
% Boundaries = bwboundaries(BW,'holes'); 
% Boundaries = Boundaries{1};
% 
% [x,y] = reducem( Boundaries(:,2), Boundaries(:,1),20);
% 
% % Find the boundary that corresponds to the Occupancy grid
% pos = SSSWindow.InteractiveLayer(x ,y);

DATA.SSS_data = DATA.SSS_data.MergeFeatures(region);

DrawSSSImage


end


%% Crop Image and Data
function CropImage(~,~)
global SSSWindow DATA

[~, rect] = SSSWindow.CropBox;

a = round(rect(1));
b = a + round(rect(3));
b = min(b, size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Image, 2));

fields = fieldnames(DATA.SSS_data.Images);

for ii=1: length(fields)
    
    DATA.SSS_data.Images(SSSWindow.mouseCounter).(fields{ii}) = DATA.SSS_data.Images(SSSWindow.mouseCounter).(fields{ii})(:,a:b);
 
end

DrawSSSImage()                                                              % Display SSS Image

end


%% Rotate Image and Data 180
function FlipImage(~,~ )
global DATA SSSWindow

fields = fieldnames(DATA.SSS_data.Images);

for ii=1: length(fields)
    
    DATA.SSS_data.Images(SSSWindow.mouseCounter).(fields{ii}) = flipud(DATA.SSS_data.Images(SSSWindow.mouseCounter).(fields{ii}));
 
end

DrawSSSImage()                                                              % Display SSS Image

end


%% Slice SSS Image into smaller parts
function DiceImage(~,~,panelName)
global SSSWindow DATA

[row,col] = size(DATA.SSS_Data(SSSWindow.mouseCounter).image);                       % Size ofthe Sonar Image

height = 256;                                                               % Divide the image along the sonar trench

if height > col                                                             % Make Sure that the image is large enough
    disp("Image is too small to divide")
    return
end

if height > 300                                                             % In case we forget to reduce the image size
    ReduceImage(0,0);
end


numImages = ceil(col/height);                                               % Number of herizontal images to be produced

compundLength = numImages * height;                                         % Total length of the sub images

overlap = floor((compundLength - col) / (numImages -1));                    % Amount of overlap between sub images

step = height - overlap - 1;                                                % Indexing step

start = 1;
stop = height;

figure('Name','Diced Sonar Images','NumberTitle','off')

x = 1:height;
[xx,yy] = meshgrid(x,x);
inter = 1;
for ii = 1: 2: numImages*2
    
    DATA.SSS_Data(SSSWindow.mouseCounter).subImage{ii}   = DATA.SSS_Data(SSSWindow.mouseCounter).image(1:height,    start:stop);
    DATA.SSS_Data(SSSWindow.mouseCounter).subImage{ii+1} = DATA.SSS_Data(SSSWindow.mouseCounter).image(row-height+1:row,start:stop);
    
    if isfield(DATA.SSS_Data,'mask')
        DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii}   = DATA.SSS_Data(SSSWindow.mouseCounter).mask(1:height,    start:stop);
        DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii+1} = DATA.SSS_Data(SSSWindow.mouseCounter).mask(row-height+1:row,start:stop);
    end
    
    start = start + step;
    stop  = start + height-1;
    
    subplot(2,numImages,inter)
    imshow(DATA.SSS_Data(SSSWindow.mouseCounter).subImage{ii})
    if isfield(DATA.SSS_Data,'mask')
        hold on
        s = scatter(xx(DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii}),yy(DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii}));
        s.MarkerEdgeAlpha = 0.3;
        hold off
    end
    
    subplot(2,numImages, numImages + inter )
    imshow(DATA.SSS_Data(SSSWindow.mouseCounter).subImage{ii+1})
    if isfield(DATA.SSS_Data,'mask')
        hold on
        s = scatter(xx(DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii+1}),yy(DATA.SSS_Data(SSSWindow.mouseCounter).subImage_mask{ii+1}));
        s.MarkerEdgeAlpha = 0.3;
        hold off
    end
    
    
    inter = inter + 1;
end

Panel = SSSWindow.Panels(panelName); 
Panel.button(10).Push.Visible = 'on';   

end


%% Display SSS Image in GUI window
function DrawSSSImage(fromPanel)
global MainWindow SSSWindow DATA 

if nargin > 0
    lon_lat = DATA.SSS_data.Path(SSSWindow.mouseCounter).lon_lat;           % Get Path
    MainWindow.Addlayer('Path', [lon_lat(:,2), lon_lat(:,1)], fromPanel);   % Add datapoints to the GUI
    
    [lat, lon] = DATA.SSS_data.Feature_latlon;                              % Get Feature locaiton in lat lon
    if ~isempty(lat)    
        MainWindow.Addlayer('Marker', [lon, lat], fromPanel);               % Add Features to the GUI
    end
end

x = 1: size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Image,2);
y = 1: size(DATA.SSS_data.Images(SSSWindow.mouseCounter).Image,1);
[xx, yy] = meshgrid(x,y);

SSSWindow.BackGround('LON',xx,'LAT',yy,'Image', DATA.SSS_data.Images(SSSWindow.mouseCounter).Image);
SSSWindow.Addlayer('Regions', DATA.SSS_data.Images(SSSWindow.mouseCounter).Mask, 'Transformations' );
SSSWindow.Header.String = DATA.SSS_data.RawData(SSSWindow.mouseCounter).log;

% figure
% imshow(DATA.SSS_data.Images(SSSWindow.mouseCounter).Image)

end


%% Mouse Callback finction
function MouseCallBack(input)
global SSSWindow;                                                           % Side Scan Sonar GUI Window

fromPanel = input{1};
panelName = input{2};

Panel = SSSWindow.Panels(panelName);                                        % Extract the button pane from its mapped container to access the structure

SSSWindow.DeleteLayer('OcupancyGrid');

% Update Imagry
DrawSSSImage(fromPanel)

end


%% Image Rectangle
function GetRectangularPatch(~, ~)
global SSSWindow ui DATA

pos = SSSWindow.Slid_Rect(210,210);

image = DATA.SSS_data.Images(SSSWindow.mouseCounter).Image( pos(2) : pos(2) + pos(4), pos(1) : pos(1) + pos(3) );

image = imresize(image, 0.5, 'method', 'nearest');

figure
imshow(image)

% Save Data
disp("\n\nChoose Loaction to Save image patch'")

file = ui.saveFile('image_patch','jpeg');                           % GUI to get file path to save file
if isempty(file)                                                            % Check if the Data Input box was cnaceled
   return
end

imwrite(image, file)                                                        % Save Image

end