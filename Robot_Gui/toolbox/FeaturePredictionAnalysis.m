

%% Feature Prediction Analysis Tool
function FeaturePredictionAnalysis(homePanel)
global FPAWindow;                                                           % Side Scan Sonar GUI Window

disp('Side Scan Sonar')

% Make SSS GUI Window ----------------------------------------------------------------------------------------------------------------------------------------
windowName = 'Feature Prediction Analysis';                                 % GUI Windo Name
panelName  = 'Tally';                                                       % Button Panel name

FPAWindow = GUIwindow(windowName);                                          % Instanciate the GUI Windo
FPAWindow.GUIpanels(panelName,'make','position',[0.8 0 0.2 1],'backgroundcolor',[0.92941  0.6902  0.12941]);  % Create putton panel for primary GUI functions
FPAWindow.Header.BackgroundColor = [0  1  1];

% Make buttons in the button panel: 
%{ Buttin type,   Button Name,   Callback function and their inputs,   Callback function for slider listeners }
FPAWindow.GUIbuttons( panelName,'add', ...
    {'pushbutton', 'Load Data', {@LoadData,homePanel, panelName}, [], 'Load Processed Sonar Images and their Binar Predictions';   %   1
     'space',      ' ',         [],                               [], '';                                                          % --2
     'pushbutton', '+ Corect',  @PluseCorrect,                    [], 'Tally one Corectly Marked Feature';                         %   3
     'pushbutton', '+ Wrong',   @PluseWrong,                      [], 'Tally one incorectly marked feature';                      %   4
     'pushbutton', '+ Mist',    @PluseMiss,                       [], 'Tally one missed feature';                      %   4
    });

FPAWindow = FPAWindow.AdjustAplahMax(0.75);
end


%% Get Data
function LoadData(~,~,fromPanel, panelName)
global UserInputs FPAWindow DATA 

FPAWindow.Header.String = 'Select Data Input';                              % Prompt to select data
FPAWindow.mouseCounter = 1;

fileType = 'jpeg';
fun = @(file) imread(file);

sonarPath = UserInputs.getDir('Sonar Images',fileType);                             % Get folder GUI
predPath  = UserInputs.getDir('Feature Preditions',fileType);                       % Get folder GUI

if isempty(sonarPath) || isempty(predPath)                                  % Check if the Input box was cnaceled
    disp('Get Directory Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

SonarImages = BatchDir(sonarPath,fileType,fun);                             % Get Data
featurePred = BatchDir(predPath, fileType,fun); 


if isempty(SonarImages) || isempty(featurePred)                             % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end

FPAWindow.AddMouseWheel( @MouseCallBack, {fromPanel, panelName}, 1, numel(SonarImages) ); 

DATA.FeaturePreditons.SonarImages = SonarImages;
DATA.FeaturePreditons.featurePred = featurePred;

DATA.FeaturePreditons.Correct   = 0;
DATA.FeaturePreditons.Incorrect = 0;
DATA.FeaturePreditons.Miss      = 0;

ShowImage(fromPanel,panelName)
end


%% Mouse Callback finction
function MouseCallBack(input)
ShowImage(input{1},input{2})                          % Update Imagry
end

%% Tally the Results
function PluseCorrect(~,~)
global DATA FPAWindow
DATA.FeaturePreditons.Correct = DATA.FeaturePreditons.Correct + 1;

indx = FPAWindow.mouseCounter;
FPAWindow.Header.String = sprintf('Sonar Image: %d \t Correct %d,\tIncorect: %d, \tMissed %d', indx, ... 
    DATA.FeaturePreditons.Correct, ... 
    DATA.FeaturePreditons.Incorrect, ... 
    DATA.FeaturePreditons.Miss);
end


function PluseWrong(~,~)
global DATA FPAWindow
DATA.FeaturePreditons.Incorrect = DATA.FeaturePreditons.Incorrect + 1; 

indx = FPAWindow.mouseCounter;
FPAWindow.Header.String = sprintf('Sonar Image: %d \t Correct %d,\tIncorect: %d, \tMissed %d', indx, ... 
    DATA.FeaturePreditons.Correct, ... 
    DATA.FeaturePreditons.Incorrect, ... 
    DATA.FeaturePreditons.Miss);
end


function PluseMiss(~,~)
global DATA FPAWindow
DATA.FeaturePreditons.Miss = DATA.FeaturePreditons.Miss + 1; 

indx = FPAWindow.mouseCounter;
FPAWindow.Header.String = sprintf('Sonar Image: %d \t Correct %d,\tIncorect: %d, \tMissed %d', indx, ... 
    DATA.FeaturePreditons.Correct, ... 
    DATA.FeaturePreditons.Incorrect, ... 
    DATA.FeaturePreditons.Miss);
end



%% Compare Results
function ShowImage(fromPanel,panelName)
global FPAWindow DATA

indx = FPAWindow.mouseCounter;

mask = DATA.FeaturePreditons.featurePred{indx};

mask = imresize( mask, 2.0, 'method', 'nearest')/255;

lowBound = 140;
highBound = 512*512/3;

s = sum(sum(mask));

h = ones(5,5)/25;
mask = imfilter(mask,h);

mask(mask > 0.35) = 1;

[~, regions, ~, ~] = bwboundaries(mask);

if s <= lowBound || s >= highBound
    mask = zeros(size(mask));
else
    
    stats = regionprops(regions);
    
    in = cat(1,stats.Area) > 300 & cat(1,stats.Area) < 5000 ;
    
    for ii = 1:numel(in)
        mask(regions == ii) = in(ii);
    end
    
end

image = imresize( DATA.FeaturePreditons.SonarImages{indx}, 2.0, 'method', 'nearest');
x = 1: size(image,2);
y = 1: size(image,1);
[xx, yy] = meshgrid(x,y);


FPAWindow.BackGround('LON',xx,'LAT',yy,'Image', image);
FPAWindow.Addlayer('OcupancyGrid', mask, panelName );
FPAWindow.Addlayer('Regions', regions, panelName );

FPAWindow.Header.String = sprintf('Sonar Image: %d \t Correct %d,\tIncorect: %d, \tMissed %d', indx, ... 
    DATA.FeaturePreditons.Correct, ... 
    DATA.FeaturePreditons.Incorrect, ... 
    DATA.FeaturePreditons.Miss);

end

