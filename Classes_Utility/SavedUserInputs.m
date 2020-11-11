% Jacob Anderson
% December 12th 2017

% This is a classe deffinition that will allow useres to save UI inputs and then call the inputs at a later time
% Inputs are saved to Documents/MATLAB/Script_Inputs/<Script_Name>
%   <Script_Name> is the name of the script calling this classdef
%___________________________________________________________________________________________________________________________________________
% To use this script place the following commands in your script:
% 1)    UserInputs = SavedUserInputs(mfilename);            -->  Instanciate the SavedUserInputs Class.
%                                                                * Run this command once at the begining of your script.
%                                                                * This create a file for saving data or acces previously saved data.
%                                                                * mfilename is a MATLAB command that identifies the script being run.
%                                                                  Do not change this command unless you know what you're doing.

% 2)    UserInputs.Preferences('Varialbe_ID') = variable;   -->  Add a varialbe to this instance of the class.
%                                                                * Run this command, followed by the next one, everytime you want to save a variable.
%                                                                * 'Variable ID' is an arbatrary string that is used to identify the variable.
%                                                                * variable is the input to be saved. Usualy a file path or name but can be anything.

% 3)    UserInputs.SavePreferences();                       -->  This command saves the class instance to the hard drive. The varialbes can then be accessed at a later time

% 4)    variable = UserInputs.Preferences('Variable_ID');   -->  Retrive the saved data
%__________________________________________________________________________________________________________________________________________

classdef SavedUserInputs
    
    properties (Access= private)
        DefaultOpenPath
        FilePath
        FullFileName
        Preferences
        newData
    end
    
    methods
        % Initialize the SaveUserInput Class. Check to see if inputs have already been saved ---------------------------------------------------------------------------
        function obj = SavedUserInputs(scriptName), warning on backtrace
            
            obj.newData = true;
            
            obj.FilePath     = fullfile( userpath,'.Script_Inputs');
            obj.FullFileName = fullfile( obj.FilePath, [scriptName,'_suggestedFilePaths.mat']);
            obj.Preferences  = containers.Map('KeyType', 'char' ,'ValueType','any');
            
            
            if exist(obj.FullFileName,'file')                               % If the variable exists, Load the saved settings
                disp('Loading Saved Inputs')
                load(obj.FullFileName);
                obj.Preferences = suggestions;
                
            elseif ~exist(obj.FilePath,'dir')                               % Create the directory for the saved inputs if it does not exist
                disp('Making directory for saved inputs')
                mkdir(obj.FilePath)
                
            else                                                            % The directoy exists, but preferences have not been saved before
                disp('No Saved Inputs')
                return
            end
            
        end
        
        
        % Function to indicate weatehr new data should be selescted
        function obj = NewData(obj,tf), obj.newData = tf; end
        
        % User Interfaces ==================================================================================================================================================
        % Get Directory ----------------------------------------------------------------------------------------------------------------------------------------------------
        function filepath = GetDir(obj, name, ext )
            
            
            if ~obj.newData && any(contains( keys(obj.Preferences), ['dir_',name] ))
                dir = obj.Preferences(['dir_',name]);
                
            else
                
                fprintf('\n\n--Select %s Folder--\n', name)
                
                ext = strcat('*.', ext, '*');
                
                try    dir = uigetdir(obj.Preferences(['dir_',name]),['Select ',name,' Folder']);
                catch, dir = uigetdir(ext,['Select ',name,' Folder']);
                end
                
                if dir == 0                                                 % Check if the Data Input box was cnaceled
                    fprintf('%s Directory Input Box Cancelled\n\n', name)
                    dir = [];
                    
                else
                    fprintf('Path to the %s File: %s \n\n', name, dir);
                    obj.Preferences(['dir_',name]) = dir;                   % Add file path to saved inputs
                    obj.SavePreferences();                                  % Save Input
                end
                
            end
            
            if nargout > 0
                filepath = dir;
            end
        end
        
        
        % Get File path ----------------------------------------------------------------------------------------------------------------------------------------------------
        function file = GetFile(obj, name, ext)
            
            
            if ~obj.newData && any(contains( keys(obj.Preferences), ['file_',name] ))
                filePath = obj.Preferences(['file_',name]);
                
            else
                
                fprintf('\n\n--Select %s File--\n', name)
                
                ext = strcat('*.', ext, '*');
                
                try
                    [fileName, filePath] = uigetfile( ext, ['Select ',name,' File'],obj.Preferences(['file_',name]) );
                catch
                    [fileName, filePath] = uigetfile( ext, ['Select ',name,' File'] );
                end
                
                if fileName == 0                                                % Check if the Data Input box was cnaceled
                    fprintf('%s File Input Box Cancelled\n\n', name)
                    filePath = [];
                    
                else
                    fprintf('Path to the %s File: \n%s%s \n\n', name, filePath, fileName);
                    filePath = fullfile( filePath, fileName );
                    obj.Preferences(['file_',name]) = filePath;                               % Add file path to saved inputs
                    obj.SavePreferences();                                      % Save Input
                end
            end
            
            if nargout > 0
                file = filePath;
            end
        end
        
        
        % Get Save To File path --------------------------------------------------------------------------------------------------------------------------------------------
        function filepath = SaveFile(obj, name, ext)
            
            if ~obj.newData && any( contains( keys(obj.Preferences), ['save_',name] ))
                file = obj.Preferences(['save_',name]);
                
            else
                
                fprintf('\n\n--Where to Save %s File--\n', name)
                
                ext = strcat('.', ext);
                
                try
                    [path_, name_, ext_] = fileparts(obj.Preferences(['save_',name]));
                    
                    C = strsplit(name_, '_');
                    C{1} = num2str(yyyymmdd(datetime('today')));
                    
                    str = fullfile(path_, [strjoin(C, '_'), ext_]);
                   
                    
                    [fileName, filePath] = uiputfile( ['*',ext], ['Where to Save ',name,' File'], str); % obj.Preferences(['save_',name])); %,num2str(yyyymmdd(datetime('today'))),'_',name,ext]);
                catch
                    [fileName, filePath] = uiputfile( ['*',ext], ['Where to Save ',name,' File'], [num2str(yyyymmdd(datetime('today'))),'_',name,ext]);
                end
                
                if fileName == 0
                    fprintf('%s Save File Input Box Cancelled \n\n', name)
                    file = [];
                    
                else
                    fprintf('%s File Saved To: \n%s%s \n\n', name, filePath, fileName);
                    file = fullfile(filePath,fileName);
                    obj.Preferences(['save_',name]) = file;                  % Add file path to saved inputs
                    obj.SavePreferences();                                      % Save Input
                    
                end
                
            end
            
            if nargout > 0
                filepath = file;
            end
            
        end
        
        
        % Get Geotiff file -------------------------------------------------------------------------------------------------------------------------------------------------
        function gt = GetGeoTiff(obj,varargin)
            
            gt = Geotiff;
            
            if ~obj.newData && any( contains( keys(obj.Preferences), 'save_geo-tiff' ))
                file = obj.Preferences('save_geo-tiff');
                gt = gt.FromFile(file);
                return
                
            elseif ~obj.newData && any( contains( keys(obj.Preferences), 'file_geo-tiff' ))
                file = obj.Preferences('file_geo-tiff');
                gt = gt.FromFile(file);
                
            else
                
                answer = questdlg('Where would you Like to get the Geo-tiff?', ...
                                  'Source Geo-tiff', ...
                                  'File on Computer','Google', 'Skip Geo-tiff', ...
                                  'File on Computer');
                
                %            answer = obj.GetTiffGUI();
                
                % Handle response
                switch answer
                    
                    case 'Skip Geo-tiff'   % User Does not want a Geotiff ---------------------------------------------------------
                        fprintf('\nSkiping Geotiff\n\n')                                     % End script if there isn't an input for it to use
                        return
                        
                        
                    case 'File on Computer'     % Get Geotiff from existing file -----------------------------------------------------------
                        
                        file = obj.GetFile('geo-tiff','tif');                   % GUI to get the name and file path of a file
                        
                        if isempty(file), return, end                           % End script if there isn't an input for it to use
                        
                        gt = gt.FromFile(file);
                        
                        return
                        
                        
                    case 'Google'   % Get new geotiff from google ----------------------------------------------------------------------------
                        
                        if nargin < 2                                           % Make sure a center coordinate has been specified for the geotiff
                            disp('A latitude and longitude position must be supplied for the center of the geotiff')
                            return
                        end
                        
                        lon = mean(varargin{1});
                        lat = mean(varargin{2});
                        
                        zoom = inputdlg('Zoom level between 1 and 18','Zoom', [1 10], {'18'});
                        zoom = str2double(zoom);
                        
                        gt = gt.FromGoogle(lat, lon, zoom);
                        
                        
                    case 'Sentinel-2 Cloudless'   % Get new geotiff from google ----------------------------------------------------------------------------
                        
                        if nargin < 2                                           % Make sure a center coordinate has been specified for the geotiff
                            disp('A latitude and longitude position must be supplied for the center of the geotiff')
                            return
                        end
                        
                        lon = [min(varargin{1}), max(varargin{1})];
                        lat = [min(varargin{2}), max(varargin{2})];
                        
                        gt = gt.FromWMS(lat, lon);
                        
                        
                    otherwise
                        fprintf('\nSkiping Geotiff\n\n')                                     % End script if there isn't an input for it to use
                        return
                        
                end
                
                fig = figure('Name','Geottiff','NumberTitle','off');
                gt.Show
                
                file = obj.SaveFile('geo-tiff','tif');                  % GUI to get file path to save file
                fprintf('\nSave File Path and name: %s \n\n', file);
                
                if isempty(file)                                        % Check if the Data Input box was cnaceled
                    disp('Input box Cancelled, Geotiff not saved!!')
                else
                    gt.ExportTiff(file)
                end
                
                close(fig)
                
            end
        end
        
        
        % Manage API Keys
        function apiKey = getAPIkey(obj, webservice )
            
            % first run, check if API key file exists ------------------------------------------------------
            apiFilePath = fullfile(userpath,'.ApiKeys', 'Keys.mat'); % Check for API Key directory
            
            if ~exist(apiFilePath,'file')
                obj.MakeKeyChain(apiFilePath);                     % Make directory if it does not exist
            end
            
            load(apiFilePath,'Keys');
            
            if isKey(Keys, webservice )
                apiKey = Keys(webservice);
                
            else
                
                disp('-----------------------------------------------------')
                disp('  Can Not Retrive Google Map Without API Key')
                disp('  API Key can be obtained for free at:')
                disp('  https://cloud.google.com/maps-platform/')
                disp("  Goto this web address and click on 'Get Started'")
                disp('-----------------------------------------------------')
                
                title = [ webservice, ' API Key Required'];
                prompt = {'Enter API Key'};
                dims = [1 60];
                definput = {'api-key'};
                answer = inputdlg(prompt,title,dims,definput);
                
                if isempty(answer)
                    apiKey = [];
                    return
                end
                
                apiKey = answer{1};
                
                Keys(webservice) = apiKey;
                
                save(apiFilePath,'Keys')
                
            end
        end
        
        
        % Settings for other types of inputs -----------------------------------------------------------------------------------------------------------------------------
        function settings = GetSettings(obj, type)
            
            if isKey(obj.Preferences, ['public',type] )
                settings = obj.Preferences(['public',type]);
            else
                settings = [];
            end
        end
        
        
        % Save settings for other types of inputs -----------------------------------------------------------------------------------------------------------------------------
        function SaveSettings(obj, type, vars)
            obj.Preferences(['public',type]) = vars;
            obj.SavePreferences;
        end
        
        
        function WaypointList(obj, fig)
            
            if nargin > 1, [x,y] = getpts(fig);
            else, [x,y] = getpts;
            end
            
            writematrix(obj.GetFile('Waypoint_List','txt'), [x,y]);
        end
        
    end
    
    
    methods (Access= private)
        
        % Save User Inputs -------------------------------------------------------------------------------------------------------------------------------------------------
        function SavePreferences(obj,varargin)
            suggestions = obj.Preferences;
            save(obj.FullFileName,'suggestions');
            
        end
        
        
        % Make Directory for API keys --------------------------------------------------------------------------------------------------------------------------------------
        function MakeKeyChain(~,filename)
            
            [filepath,~,~] = fileparts(filename);
            
            fprintf('\nCreating diretory for API Keys:%s \n\n', filepath)
            
            mkdir(filepath);
            
            Keys = containers.Map('KeyType', 'char' ,'ValueType','char');
            
            save(filename,'Keys');
            
        end
    
        
        % Set Default Open Path -------------------------------------------------------------------------------------------------------------------------------------------------
        function obj = SetDefaultOpenPath(obj,~)
            
            parts = strsplit(pwd,{'/','\'});
            
            if ispc
                userHome = ['\',parts{2}];
            else
                userHome = ['/',parts{2}];
            end
            
            filelist=dir(fullfile(userHome, '**', '_WrokingData'));
            
            if isempty(filelist) || any( strcmpi(parts, 'trash'))
                obj.DefaultOpenPath = pwd;
            else
                obj.DefaultOpenPath = filelist(1).folder;
            end
            
        end
        
        
        % Set Default Save Path -------------------------------------------------------------------------------------------------------------------------------------------------
        function obj = SetDefaultSavePath(obj,~)
            
            parts = strsplit(pwd,{'/','\'});
            
            if ispc
                userHome = ['\',parts{2}];
            else
                userHome = ['/',parts{2}];
            end
            
            filelist=dir(fullfile(userHome, '**', '_OutPuts'));
            
            if isempty(filelist) || any( strcmpi(parts, 'trash'))
                obj.DefaultSavePath = pwd;
            else
                obj.DefaultSavePath = filelist(1).folder;
            end
            
        end
        
        
        % Get Geotiff GUI -------------------------------------------------------------------------------------------------------------------------------------------------
        function answer = GetTiffGUI(obj)
            
            fig = figure('Name',figName, 'units','normalized',...
                'Position',[0.1  0.1   0.7  0.7],...
                'Color',[0.8  0.8  0.8],...
                'NumberTitle','off',...
                'MenuBar','none',...
                'ToolBar','none' );
            
            c = uicontextmenu;
            uimenu(c,'Label','Save Figure','Callback',@obj.MakeFigure);
            
            
            axes(fig, ...
                'Position',[0.1  0.1   0.65 0.7], ...
                'Color',[0.25, 0.25, 0.25], ...
                'UIContextMenu', c);
            
            uicontrol(fig,'Style','text','units','normalized', ...
                'Position',[0  0.95  0.8  0.05], ...
                'BackgroundColor',[0.3020  0.7490  0.9294], ...
                'String', 'Hello', ...
                'FontName', 'FixedWidth', ...
                'FontSize', 14, ...
                'FontWeight','bold', ...
                'HorizontalAlignment', 'left' );
            
            answer = [];
        end
        
        
        function obj = MakeFigure(obj)
            
            % Get Geotiff GUI goes here
        end
        
    end
    
end



