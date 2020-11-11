% Side Scan Sonar Data Class
% Jacob Anderson
% 8/28/2019

classdef SSS_Data
    
    properties
        RawData
        Images
        Features
        Path
        EM_Index
        GlobalMap
        SubMap
    end
    
    properties (Access = private)
        EM_Alt
        SonarRanges
        bathymetry
        portImages
        starImages
        max_range
    end
    
    
    
    %% General Methods ----------------------------------------------------
    methods
        
        % Instanciate Class by loading Raw Data
        function obj = SSS_Data(file), warning on backtrace
            
            % Do stuff with matfile()s
            
            if nargin < 1
                obj.RawData       = struct;
                obj.Images.Image  = [];
                obj.Images.Mask   = [];
                obj.Images.Trough = [];
                obj.Images.XX     = [];
                obj.Images.YY     = [];
                obj.Path.utm      = [];
                obj.Path.lon_lat  = [];
                
            else
                if isfile(file)                                             % Field is an individual file
                    
                    obj.RawData = Logdoc2mat(file);                         % Load Data
                    
                elseif isfolder(file)                                       % Feild is a directory
                    
                    fileType = 'logdoc';
                    fun = @(file) Logdoc2mat(file);
                    
                    obj.RawData = BatchDir(file,fileType,fun);              % Load Data
                    
                else
                    disp("EM_Data --> Invalid File Path")                   % Error State
                    obj = [];
                    return
                end
                
                if isempty(obj.RawData)
                    disp("Empty Sonar Data Structure")
                    disp("Doulbe check the directroy / file that you chose")
                    obj = [];
                end
                
                % Fill in Sonar Track Data
                for ii = 1: size(obj.RawData)
                    
                    % Create Sonar images for each track
                    portSonar = cat(1,obj.RawData(ii).portSonar);           % Port Sonar data structure
                    starSonar = cat(1,obj.RawData(ii).starSonar);           % Starboard Sonar data structure
                    
                    portEcho = cat(2,portSonar.EchoStrength);               % Concatinate the image date from the sonar data structures
                    starEcho = cat(2,starSonar.EchoStrength);
                    
                    obj.Images(ii).Image  = [flipud(portEcho); starEcho];   % Create FullSonar Image
                    
                    % Create Corespnomding image mask, sonar trough indicator and geo-mesh
                    obj.Images(ii).Mask   = zeros(size(obj.Images(ii).Image));
                    obj.Images(ii).Trough = false(size(obj.Images(ii).Image));
                    obj.Images(ii).XX     = zeros(size(obj.Images(ii).Image));
                    obj.Images(ii).YY     = zeros(size(obj.Images(ii).Image));
                    
                    % Path Data in utm and lat-lon
                    utm     = obj.RawData(ii).position;                     % Vehicle location in UTMs [ Easting, Northing, Altitude]
                    utmzone = obj.RawData(ii).utmZone;                      % UTM-zone: Longitude and Latitude bands
                    
                    [lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);
                    
                    obj.Path(ii).utm     = utm;
                    obj.Path(ii).utmZone = utmzone;
                    obj.Path(ii).lon_lat = [lon, lat];
                    
                end
            end
            
            obj.EM_Index = [];
            obj.EM_Alt   = [];
            obj.GlobalMap = [];
            
        end
        
        
        % Time Index The Side Scan Sonar Data to Ecomapper Data
        function obj = Index2EM_Data(obj, em_times)
            
            sizes = cellfun(@(x)size(x,1),{obj.RawData.timeStamp});         % Get Max Size of data fileds
            sizes = max(sizes);
            
            obj.EM_Index = zeros(sizes, size(obj.RawData,1));               % Pre-alocate memory
            
            for ii = 1:size(obj.RawData,1)                                  % Iterate through number of fields
                
                sss_times = obj.RawData(ii).timeStamp;
                
                for jj = 1: size(sss_times,1)                               % Iterate through entires in each field
                    
                    diffs = abs(sss_times(jj) - em_times);
                    ind = find( diffs == min(diffs(:)));
                    obj.EM_Index(jj,ii) = ind(1);
                    
                end
                
            end
        end
        
        
        % Get Coresponding Altimiter Readings
        function obj = Add_EM_Alt(obj, em_alt)
            
            [r,c] = size(obj.EM_Index);                                     % Get size of indexing matrix
            
            obj.EM_Alt = zeros(r,c);                                        % Pre-alocate memory
            
            for ii = 1:c                                                    % Iterate through number of fields
                for jj = 1:r                                                % Iterate through entries in each field
                    
                    ind = obj.EM_Index(jj,ii);
                    if ind == 0                                             % Zero means there is no mor data
                        break
                    end
                    
                    obj.EM_Alt(jj,ii) = em_alt(ind);                        % Assign altimiter data
                    
                end
            end
            
            obj = obj.GeoMesh_From_Altimeter;                                        % Create Goe-Mesh using altimiter readings
        end
        
        
        % Establish Sonar Ranges
        function obj = MakeSoanrRanges(obj,maxRange)
            
            height = size(obj.Images(1).Image,1);                           % Dimentions os the sonar track
            midLine = height/2;
            s_dist = 1:height;                                              % Sonar Ranges [pixels]
            obj.SonarRanges = (s_dist - midLine) .* maxRange/midLine;       % Sonar Ranges [meters]
            obj.max_range = maxRange;
        end
        
        
        % Get the Sonar Scanline
        function [sx, sy, z] = ScanLine(obj, xy, theta, depth)
            
            theta = theta + 90;                                             % Orientation of the scan line [deg]
            
            theta( theta < 0 )   = theta( theta < 0 ) + 360;                % Keep  360 >= theta >= 0
            theta( theta > 360 ) = theta( theta > 360 ) - 360;
            
            if nargin < 4
                sx = xy(:,1) + obj.SonarRanges .* sind( theta );        % Posible UTMs of the scan line
                sy = xy(:,2) + obj.SonarRanges .* cosd( theta );
                return
                
            elseif isa(depth,'Bathymetry_Map')
                
                px = xy(:,1) + obj.SonarRanges .* sind( theta );        % Posible UTMs of the scan line
                py = xy(:,2) + obj.SonarRanges .* cosd( theta );
                
                depths = depth.Depth([px',py']);
                
                hypot = sqrt( obj.SonarRanges'.* obj.SonarRanges' + depths.*depths) .* sign(obj.SonarRanges');
                
                hypot_ = hypot;
                hypot_(hypot_ == 0) = [];
                min_ = min(abs(hypot_));
                
                sx = NaN(size(px));
                sy = NaN(size(py));
                sz = NaN(size(py));
                
                for a = 1: numel(obj.SonarRanges)
                    
                    if abs(obj.SonarRanges(a)) < min_ , continue, end       % Skip sonar Trought
                    
                    diffs = abs(obj.SonarRanges(a) - hypot);
                    ind = find( diffs == min(diffs(:)));
                    
                    sx(a) = px(ind(1));                                     % Global position of Sonar Pixels
                    sy(a) = py(ind(1));
                    sz(a) = depths(ind(1));
                end
                
                if nargout > 2
                    z = sz;
                end
                
            else
                dd = obj.SonarRanges.*obj.SonarRanges - depth*depth;
                dx = sqrt(abs(dd)) .* sign(obj.SonarRanges);
                
                sx = xy(:,1) + dx .* sind( theta(jj) );                     % Posible UTMs of the scan line
                sy = xy(:,2) + dx .* cosd( theta(jj) );
                
                sx(dd <=0 ) = NaN;
                sy(dd <=0 ) = NaN;
            end
            
            
        end
        
        
        % Create Geomesh from Altimiter
        function obj = GeoMesh_From_Altimeter(obj)
            
            for ii = 1: size(obj.RawData)                                   % Iterate through th evarious sonar tracks
                
                [height,width] = size(obj.Images(ii).Image);                % Dimentions os the sonar track
                
                utm     = obj.RawData(ii).position(:,1:2);                  % Get Data from the sonar track
                heading = obj.RawData(ii).heading;
                depth   = obj.EM_Alt(1:width ,ii) * 1.2;
                
                midLine = height/2;
                
                s_dist = 1:height;                                          % Sonar Ranges [pixels]
                s_dist = (s_dist - midLine) .* 30/midLine;                  % Sonar Ranges [meters]
                
                for jj = 1:width
                    
                    h_dist = s_dist.*s_dist - depth(jj)*depth(jj);          % Horizontal distance along the seafloor, radiating out from the vehicle [meters]
                    
                    trought = h_dist <= 0;                                  % -number indicate the sonar trought
                    
                    h_dist(trought) = 0;                                    % Clear Sonar Trough
                    h_dist = sqrt(h_dist) .* sign(s_dist);                  % Horizontal distance along the seafloor [meters]
                    
                    xy = utm(jj,:);                                         % Location of the vehicle [utms "meters"]
                    
                    theta = heading(jj) + 90;                               % Orientation of the feature [deg]
                    
                    theta( theta < 0 )   = theta( theta < 0 ) + 360;        % Keep  360 >= theta >= 0
                    theta( theta > 360 ) = theta( theta > 360 ) - 360;
                    
                    obj.Images(ii).XX(:,jj) = xy(:,1) + h_dist .* sind( theta ); % Global position of Sonar Pixels
                    obj.Images(ii).YY(:,jj) = xy(:,2) + h_dist .* cosd( theta );
                    
                    obj.Images(ii).XX(trought,jj) = NaN;                    % Global position of Sonar Pixels
                    obj.Images(ii).YY(trought,jj) = NaN;
                    obj.Images(ii).Trough(trought, jj) = true;
                    
                end
                
            end
            
        end
        
        
        % Use Bathymetry Map to Creat Geomesh
        function obj = GeoMesh_From_Bathymetry(obj, bathyMap)
            
            for ii = 1: size(obj.RawData)                                   % Iterate through th evarious sonar tracks
                
                width = size(obj.Images(ii).Image,2);                       % Dimentions os the sonar track
                
                utm     = obj.RawData(ii).position(:,1:2);                  % Get Data from the sonar track
                heading = obj.RawData(ii).heading;
                
                for jj = 1:width
                    
                    [px, py] = obj.ScanLine( utm(jj,:), heading(jj), bathyMap);
                    
                    obj.Images(ii).XX(:,jj) = px;
                    obj.Images(ii).YY(:,jj) = py;
                    
                    obj.Images(ii).Trough(isnan(obj.Images(ii).XX)) = true;
                    
                end
                
            end
            
        end
        
        
        % Feature Position in lat - lon
        function [lat, lon] = Feature_latlon(obj)
            
            features = cat(1,obj.Features.Centroid);
            
            if ~isempty(features)
                utmZone = obj.Features.UTMzone;
                [lat, lon] = utm2deg(features(:,1), features(:,2), repmat(utmZone, size(features,1),1));
            else
                lat = [];
                lon = [];
            end
            
        end
        
        % Save Features
        function SaveFeatures(obj,filePath)
            features = obj.Features;
            disp("Saving Features")
            save(filePath,'features')
            disp("  Done")
        end
        
        
        % Load Features
        function obj = LoadFeatures(obj,filePath)
            
            disp("Loading Features")
            new_features = load(filePath);
            obj.Features = new_features.features;
            
            obj = obj.Propigate_Features;
            disp("  Done")
        end
        
        
        % Save Object Function --> used by save()
        function s = saveobj(obj)
            s.Features = obj.Features;
            s.Images   = obj.Images;
            s.RawData  = obj.RawData;
            s.EM_Index = obj.EM_Index;
            s.EM_Alt   = obj.EM_Alt;
            s.Path     = obj.Path;
        end
        
    end
    
    
    methods(Static)
        
        % Load object function --> used by load()
        function obj = loadobj(s)
            if isstruct(s)
                newObj =  SSS_Data;
                newObj.Features = s.Features;
                newObj.Images   = s.Images;
                newObj.RawData  = s.RawData;
                newObj.EM_Index = s.EM_Index;
                newObj.EM_Alt   = s.EM_Alt;
                newObj.Path     = s.Path;
                
                obj = newObj;
            else
                obj = s;
            end
        end
        
    end
    
    
    
    %% ISAM Methods -------------------------------------------------------
    methods
        
        
        % --- Preprocess Sonar Data ---------------------------------------
        % Dice the entire sonar track
        function obj = DiceSonarTrack(obj, windowSize, overlap, savefile)
            
            [filepath,name,~] = fileparts( savefile );
            
            portSonar = cat(1, obj.RawData.portSonar);
            starSonar = cat(1, obj.RawData.starSonar);
            
            masterImage = [flipud(cat(2,portSonar.EchoStrength)); cat(2,starSonar.EchoStrength)]; % Put the port and starboard images together
            
            step    = windowSize * (1-overlap);                             % Distance between each image
            
            [row,col] = size(masterImage);                                  % Size ofthe Sonar Image
            
            if windowSize > col                                                 % Make Sure that the image is large enough
                disp("Image is too small to divide")
                return
            end
            
            numImages = floor(col/step);                                    % Number of herizontal images to be produced
            
            % Pre alocate memory for image pathces
            obj.portImages(numImages).image(windowSize, windowSize) = 0;
            obj.portImages(numImages).binary(windowSize, windowSize) = 0;
            obj.portImages(numImages).binName = " Names ";
            obj.portImages(numImages).imName = "Name";
            obj.portImages(numImages).index = [0,0,0,0];
            
            obj.starImages(numImages).image(windowSize, windowSize) = 0;
            obj.starImages(numImages).binary(windowSize, windowSize) = 0;
            obj.starImages(numImages).binName = " Names ";
            obj.starImages(numImages).imName = "Name";
            obj.starImages(numImages).index = [0,0,0,0];
            
            iter = 0;
            for start = 1: step: (col - windowSize)
                
                stop  = start + windowSize-1;
                iter = iter + 1;
                
                portImage = masterImage(1:windowSize,         start:stop);
                starImage = masterImage(row-windowSize+1:row, start:stop);
                
                obj.portImages(iter).image = portImage;
                obj.portImages(iter).imName = sprintf('%s_port_%03d.jpeg', name, iter);
                obj.portImages(iter).index = [1, windowSize, start, stop];
                
                obj.starImages(iter).image = starImage;
                obj.starImages(iter).imName = sprintf('%s_star_%03d.jpeg', name, iter);
                obj.starImages(iter).index = [row-windowSize+1, row, start, stop];
                
                if nargin > 3                                               % Down size and Save Image as a jpeg for Feature extraction CNN
                    starImage_ = imresize( starImage, 0.50, 'method', 'nearest');
                    portImage_ = imresize( portImage, 0.50, 'method', 'nearest');
                    
                    imwrite(starImage_, sprintf('%s/%s_star_%03d.jpeg', filepath, name, iter))     % Save Image
                    imwrite(portImage_, sprintf('%s/%s_port_%03d.jpeg', filepath, name, iter))     % Save Image
                end
            end
            
            obj.portImages(iter+1 : end) = [];                                              % Get Rid of empyt entries
            obj.starImages(iter+1 : end) = [];
            
        end
        
        
        % Get binaries for the sonar track
        function obj = GetDicedBinaries(obj,filePath)
            
            files = fullfile(filePath, 'binary_*.jpeg' );
            logfiles = dir(fullfile(files));
            fprintf(' %d binary_*.jpeg files\n',size(logfiles,1));
            
            portiter = 1;
            stariter = 1;
            
            for idy = 1: size(logfiles,1)
                
                img = single(imread(  fullfile( logfiles(idy).folder, logfiles(idy).name )  ));
                img(img < 255) = 0;
                img = img/255;
                img = imresize( img, 2.0, 'method', 'nearest');
                
                % Check that the files have data. i.e.. more than just a header
                if contains(logfiles(idy).name, "port" )
                    obj.portImages(portiter).binary = img;
                    obj.portImages(portiter).binName = logfiles(idy).name;
                    portiter = portiter + 1;
                    
                else
                    obj.starImages(stariter).binary = img;
                    obj.starImages(stariter).binName = logfiles(idy).name;
                    stariter = stariter + 1;
                    
                end
                
            end
            
        end
        
        
        % Get Sonar Features from binaries
        function obj = GetFeaturesFromBinaries(obj, maxSize, minSize, imageSize, imagePlot)
            
            portSonar = cat(1, obj.RawData.portSonar);
            starSonar = cat(1, obj.RawData.starSonar);
            
            masterImage = [flipud(cat(2,portSonar.EchoStrength)); cat(2,starSonar.EchoStrength)]; % Put the port and starboard images together
            
            masterBinary = zeros(size(masterImage));
            
            % Merge Binaries
            for ii = 1:1:size(obj.portImages,2)
                
                ind = obj.portImages(ii).index;
                masterBinary(ind(1):ind(2), ind(3):ind(4)) = masterBinary(ind(1):ind(2), ind(3):ind(4)) | obj.portImages(ii).binary;
                
                ind = obj.starImages(ii).index;
                masterBinary(ind(1):ind(2), ind(3):ind(4)) = masterBinary(ind(1):ind(2), ind(3):ind(4)) | obj.starImages(ii).binary;
                
                
            end
            
            filter = ones(5,5)/25;
            
            for ii = 1:10
                masterBinary = imfilter(masterBinary,filter);
                masterBinary(masterBinary > 0.35) = 1;
            end
            
            %             trought = cat(2, obj.Images.Trough);
            %             masterBinary(trought) = 0;
            
            [boundaries, regions, ~, ~] = bwboundaries(masterBinary);
            
            stats = regionprops(regions);
            
            in =  minSize < cat(1,stats.Area) & cat(1,stats.Area) < maxSize;
            
            boundaries = boundaries(in);
            stats = stats(in);
            
            q = numel(boundaries);
            
            % Create Structure for Sonar Feature data
            obj.Features(q).ImagePatch(imageSize, imageSize) = 0;
            obj.Features(q).ImageName = 'Name Goes Here';
            obj.Features(q).Centroid = [0,0];
            obj.Features(q).Boundary = [];
            obj.Features(q).Cov      = zeros(2,2);
            obj.Features(q).Radius   = 0;
            obj.Features(q).Index     = [0,0];
            obj.Features(q).Timestamp = datetime;
            obj.Features(1).Probs = [];
            obj.Features(1).Truth = false;
            obj.Features(1).Header = [ "ImagePatch: An image patch centered on the feature";
                " ImageName: Name of the saved image patch jpeg";
                "  Centroid: UTMs of the feature's locaiton";
                "  Boundary: Ploygon that outlines the feature in the master Image";
                "       Cov: Covariance of the Feature's locaiton relative to the vehicle";
                "    Radius: Feature's distance away from the vehicle; + is satboard, - is port";
                "     Index: Center of the feature on the master image";
                " TimeStamp: Time stamp of the feature's Index";
                "     Prods: Feature matching probabilites";
                "     Truth: Indicates ture feature matches"];
            
            sw = SlidingWindow(imageSize, size(masterImage));
            
            XX = cat(2,obj.Images.XX);
            YY = cat(2,obj.Images.YY);
            
            timeStamp = cat(1, obj.RawData.timeStamp);
            
            for ii = 1: q
                
                center = fliplr( round( stats(ii).Centroid ) );               % Centroid --> row, col  (After flip)
                [row,col] = sw.D2(center);
                
                obj.Features(ii).Centroid   = [XX(center(1),center(2)), YY(center(1),center(2))];                             % Create Structure for Sonar Feature data
                obj.Features(ii).Boundary   = boundaries{ii};
                obj.Features(ii).Cov        = zeros(2,2);
                obj.Features(ii).Radius     = obj.SonarRanges(center(1));
                obj.Features(ii).ImagePatch = masterImage(row,col);
                obj.Features(ii).Index      = center;
                obj.Features(ii).Timestamp  = timeStamp(center(2));
                obj.Features(ii).ImageName  = sprintf('Feature_%2d', ii);
                
                if nargin > 4 && imagePlot
                    figure
                    imshow(masterImage(row,col));
                end
                
            end
            
        end
        
        
        % Save Sonar Feature Images
        function obj = ExportFeatureImages(obj, savePath)
            
            [path,name,~] = fileparts(savePath);
            
            for ii = 1: size(obj.Features,2)
                
                image = imresize(obj.Features(ii).ImagePatch, 0.5, 'method', 'nearest');
                imwrite(image, fullfile(path, sprintf("%s_instance_%02d.jpeg",name,ii)))
                obj.Features(ii).ImageName = sprintf("%s_instance_%02d.jpeg",name,ii);
                
            end
        end
        
        
        % Get Feature Matching Data
        function obj = GetFeatureMatchingData(obj, file)
            
            fileID = fopen(file,'r');
            matchData = textscan(fileID, '%s %s %f', 'Delimiter',',');
            fclose(fileID);
            
            names1 = matchData{1};
            names2 = matchData{2};
            prob   = matchData{3};
            
            %             names = unique([matchData{1}; matchData{2}]);
            names = cat(1,obj.Features.ImageName);
            
            Prob = zeros(size(names,1));
            
            for ii = 1 : size(matchData{1})
                
                row = find(strcmp(names1(ii), names));
                col = find(strcmp(names2(ii), names));
                
                Prob(row, col) = prob(ii);
                Prob(col, row) = prob(ii);
                
            end
            
            obj.Features(1).Probs = Prob;
            obj.Features(1).Truth = false(size(Prob));
            
        end
        
        
        % Get EM timstamp index of features
        function obj = EM_Feature_Index(obj)
            
            emIndx = cat(1, obj.EM_Index);
            emIndx( emIndx == 0 ) = [];
            
            for ii = 1: size(obj.Features,2)
                obj.Features(ii).EM_index = emIndx(obj.Features(ii).Index(2));
            end
        end
        
        
        % Add Bathymetry Map
        function obj = Add_Batymetry(obj, bathy)
            obj.bathymetry = bathy;   
        end
        
        %--- Sonar Feature Function ---------------------------------------
        % Feature Covariance
        function obj = Feature_Make_Cov(obj, ind, heading, alt, v_cov)
            
            in = ismember([obj.Features.EM_index], ind);
            
            obj = obj.FeatureCovariance(in, heading, alt);
            
            if nargin < 5, return, end
            
            ind = find(in);
            for ii = 1:sum(in)
                obj.Features(ind(ii)).Error = obj.Features(ind(ii)).Cov + v_cov;
            end
        end
        
              
        % Return SLAM Feature
        function feature = ISAM_feature(obj, ind, xy, cov)
            
            in = ismember([obj.Features.EM_index], ind);
            
            feature(sum(in)) = sssFeature;
            
            ind = find(in);
            for ii = 1: sum(in)
                feature(ii).name        = obj.Features(ind(ii)).ImageName;
                feature(ii).xyz         = []; %[obj.Features(ind(ii)).Centroid, 0];
                feature(ii).radius      = obj.Features(ind(ii)).Radius;
                feature(ii).cov         = obj.Features(ind(ii)).Cov;
                feature(ii).times       = obj.Features(ind(ii)).Timestamp;
                feature(ii).error       = obj.Features(ind(ii)).Error;
                feature(ii).image       = obj.Features(ind(ii)).ImagePatch;
                feature(ii).vehicle     = xy;
                feature(ii).vehicle_cov = cov;
                feature(ii).ind         = ind;
            end
            
            feature(1).matchXY    = [];
            feature(1).matchName  = '';
            feature(1).matchTime  = datetime;
            feature(1).matchImage = [];
        end
        
        
        % Feature Centroid from Bathymetry map
        function feature = Feature_CentroidfromBathy(obj, feature, xy, theta)
            
            [sx, sy, z] = obj.ScanLine(xy, theta, obj.bathymetry);
            
            for ii = 1: max(size(feature))
                diffs = abs(obj.SonarRanges - feature(ii).radius);
                in = find( diffs == min(diffs(:)));
                
                feature(ii).xyz = [sx(in(1)), sy(in(1)), z(in(1))];
            end
        end
         
        % Get Sub-Map Indecies
        function obj = GetSubMapIndex(obj)
            
            obj.SubMap.Indx = [0, max(obj.EM_Index,[],1)];
            
            obj.SubMap.features = struct('name', '',       'xyz', [],         'radius', [],    ...
                                         'cov', [],        'times', datetime, 'error', [],     ...
                                         'matchXY', [],    'matchName', '',   'matchTime', [], ...
                                         'matchImage', [], 'image', []);
        end
        
    end
    
    
    methods (Access = private)
        
        % Get matching Probability ----------------------------------------------------------------------
        function prob = GetMatchingProb(obj, name1, name2 )
            
            if isempty(name1) || isempty(name2)
                prob = 0;
                return
            end
            
            prob = zeros(size(name1,1),size(name2,1));
            
            for ii = 1: size(name1,1)
                for jj= 1: size(name2,1)
                    
                    try
                    [row, col] = obj.GetProbsIndex(name1(ii,:), name2(jj));
                    prob(ii,jj) = obj.Features(1).Probs(row,col);
                    catch
                        disp("problem")
                    end
                end
            end
            
        end
        
        
        % Get Ground Truth Matching data ----------------------------------------------------------------
        function truth = GetMatchingTruth(obj, name1, name2 )
            
            if isempty(name1) || isempty(name2)
                truth = false;
                return
            end
            
            truth = false(size(name1,1),size(name2,1));
            
            for ii = 1: size(name1,1)
                for jj= 1: size(name2,1)
                    
                    [row, col] = obj.GetProbsIndex(name1(ii,:), name2(jj,:));
                    truth(ii,jj) = obj.Features(1).Truth(row,col);
                    
                end
            end
            
        end
        
        
        % Get index of probabilty matrix ------------------------------------------------------------------
        function [row, col] = GetProbsIndex(obj, name1, name2)
            
            names = cat(1,obj.Features.ImageName);
            
            row = find(strcmp(name1, names));
            col = find(strcmp(name2, names));
        end
        
        
        % Calculate Feature Covariance
        function obj = FeatureCovariance(obj, ind, theta, depth)
            
            radius = [obj.Features(ind).Radius];
            
            theta = theta + 90;                                             % Orientation of the feature [deg]
            
            theta( theta < 0 )   = theta( theta < 0 ) + 360;                % Keep  360 >= theta >= 0
            theta( theta > 360 ) = theta( theta > 360 ) - 360;
            
            % Covariance of Features
            Xerr = (abs(radius) .* sind(80) - depth .* sind(20)) .* sign(radius);
            
            mu    = arrayfun( @(x)   [x,0],       Xerr/2,                                   'UniformOutput',false);
            sigma = arrayfun( @(x)   [x,0;0,0.1], abs(Xerr/6),                              'UniformOutput',false);
            Rz    = arrayfun( @(yaw) [sind(yaw)  -cosd(yaw);  cosd(yaw)  sind(yaw)], theta, 'UniformOutput',false);
            
            xy_   = cellfun( @(x,y) mvnrnd(x,y, 100), mu, sigma, 'UniformOutput',false);
            
            Rz = repmat(Rz,(1:max(size(xy_))));
            
            xy_   = cellfun( @(x,r) RotZ(x,r), xy_, Rz,       'UniformOutput', false);
            f_cov = cellfun( @(x)   cov(x(:,1), x(:,2)), xy_, 'UniformOutput', false);
            
            if any(isnan([f_cov{:}]))
                disp("NaNs!!!! --> Paused for human inspection")
                pause
            end
            
            in = find(ind);
            for ii = 1: sum(ind) , obj.Features(in(ii)).Cov = f_cov{ii}; end
            
        end
    end
    
    
    %% GUI Methods --------------------------------------------------------
    methods
        
        % Reset
        function obj = ResetData(obj)
            for ii = 1: size(obj.Images,2)
                obj.Images(ii).Mask  = zeros(size(obj.Images(ii).Image));
            end
            
            for ii =1: size(obj.Features,2)
                obj.Features(ii).Boundary   = [];
                obj.Features(ii).Centroid   = [];
                obj.Features(ii).ImagePatch = {};
                obj.Features(ii).Sheet      = {};
            end
        end
        
        
        % Propigate Feature Masks
        function obj = AddFeatue(obj, pos, sheet)
            
            % Get Number of Features
            if isempty(obj.Features(1).Boundary)
                ii = 1;
            else
                ii = size(obj.Features,2) + 1;
            end
            
            % Creat Mask of the area inside the polygon
            [m,n] =size(obj.Images(sheet).Image);
            pos = round(pos);
            
            x = 1:n;
            y = 1:m;
            
            [XX, YY] = meshgrid(x,y);
            
            in = inpolygon(XX, YY, pos(:,1), pos(:,2));
            in = in  & ~obj.Images(sheet).Trough;
            
            obj.Images(sheet).Mask(in) = ii;
            
            % Get Centroid and boundary
            ind = sub2ind([m,n], pos(:,2), pos(:,1));
            
            XX = obj.Images(sheet).XX(ind);
            YY = obj.Images(sheet).YY(ind);
            
            % Get Image Patch
            c = round(mean(pos(:,1)));
            r = round(mean(pos(:,2)));
            
            window = SlidingWindow(210, [m,n]);
            
            [rows, cols] = window.D2([r,c]);
            
            image = obj.Images(sheet).Image(rows, cols);
            
            %             figure
            %             imshow(image)
            
            obj.Features(ii).Boundary = [XX , YY];
            obj.Features(ii).Centroid = [nanmean(XX), nanmean(YY)];
            obj.Features(ii).ImagePatch{1} = image;
            obj.Features(ii).Sheet{1} = obj.RawData(sheet).log;
            
            % Propigate the feature into all of the tracks it is present on
            obj = obj.Propigate_Features(ii, sheet);
            
        end
        
        
        % Adjust Feature
        function obj = AdjustFeatue(obj, pos, feature, sheet)
            
            % Remove Old location of feature from mask
            mask = obj.Images(sheet).Mask;
            obj.Images(sheet).Mask(mask == feature) = 0;
            
            % Creat Mask of the area inside the polygon
            [m,n] =size(obj.Images(sheet).Image);
            pos = round(pos);
            
            x = 1:n;
            y = 1:m;
            
            [XX, YY] = meshgrid(x,y);
            
            in = inpolygon(XX, YY, pos(:,1), pos(:,2));
            in = in  & ~obj.Images(sheet).Trough;
            
            obj.Images(sheet).Mask(in) = feature;
            
            % Get Centroid and boundary
            ind = sub2ind([m,n], pos(:,2), pos(:,1));
            
            XX = obj.Images(sheet).XX(ind);
            YY = obj.Images(sheet).YY(ind);
            
            % Get Image Patch
            c = round(mean(pos(:,1)));
            r = round(mean(pos(:,2)));
            
            window = SlidingWindow(210, [m,n]);
            
            [rows, cols] = window.D2([r,c]);
            
            image = obj.Images(sheet).Image(rows, cols);
            
            %             figure
            %             imshow(image)
            feature_Boundary = [obj.Features(feature).Boundary; XX, YY];
            
            feature_Boundary( isnan(feature_Boundary(:,1)), :) = [];         % Get Rid of NaNs
            feature_Boundary( isnan(feature_Boundary(:,2)), :) = [];
            
            k = boundary(feature_Boundary(:,1), feature_Boundary(:,2));
            
            obj.Features(feature).Boundary = feature_Boundary(k,:);
            
            % Average the  old and new centers
            obj.Features(feature).Centroid  = [nanmean(feature_Boundary(:,1)), nanmean(feature_Boundary(:,2))];
            
            indx = ismember(obj.Features(feature).Sheet, obj.RawData(sheet).log);
            
            obj.Features(feature).ImagePatch{indx} = image;
            
            
        end
        
        
        % Merge Features together
        function obj = MergeFeatures(obj, features)
            
            % Calculate Centroid
            feature_centroid = cat(1, obj.Features(features).Centroid);
            feature_centroid =[ mean(feature_centroid(:,1)), mean(feature_centroid(:,2)) ];
            
            % Reform the boundary
            feature_Boundary = cat(1, obj.Features(features).Boundary);
            k = boundary(feature_Boundary(:,1), feature_Boundary(:,2));
            feature_Boundary = feature_Boundary(k,:);
            
            % Merge Image Patches and sheets
            imagePatches = [obj.Features(features).ImagePatch];
            imageSheets  = [obj.Features(features).Sheet];
            
            [imageSheets, ia, ~] = unique(imageSheets);
            
            imagePatches = imagePatches(ia);
            
            % Delet the old features
            obj.Features(features) = [];
            
            % Add mergerd Feature
            ind = size(obj.Features,2) + 1;
            obj.Features(ind).Centroid = feature_centroid;
            obj.Features(ind).Boundary = feature_Boundary;
            obj.Features(ind).ImagePatch = imagePatches;
            obj.Features(ind).Sheet = imageSheets;
            
            % Propigate the changes through the data set
            obj = obj.ResetMasks;
            obj = obj.Propigate_Features;
            
        end
        
        
        % Remove Feature from Curent sheet
        function obj = RemoveFeature(obj, feature, sheet)
            
            % Remove Old location of feature from mask
            mask = obj.Images(sheet).Mask;
            obj.Images(sheet).Mask(mask == feature) = 0;
            
            % Remove Image Patch
            indx = ismember(obj.Features(feature).Sheet, obj.RawData(sheet).log);
            obj.Features(feature).ImagePatch(indx) = [];
            obj.Features(feature).Sheet(indx) = [];
            
        end
        
        
        % Get a list of the files contained in the features
        function list = GetFeatureSheetList(obj)
            
            list = {};
            
            for ii = 1:size(obj.Features,2)
                sheets = cat(1,obj.Features(ii).Sheet)';
                list = [list;sheets];
            end
            
            list = cellfun( @(x) strsplit(x, "_WP"), list, 'UniformOutput', false );
            list = cellfun( @(x) x(1), list, 'UniformOutput', false );
            list = unique([list{:}]);
            
        end
        
        
        % Save the image Patches into a file directory
        function ExportImagePatches(obj, dir, list)
            
            num_images = arrayfun( @(x) size(x.ImagePatch,2), obj.Features );
            
            count = 1;
            for ii = 1: size(obj.Features,2)
                
                if size(obj.Features(ii).ImagePatch) ~= size(obj.Features(ii).Sheet)
                    disp("Miss Matched Image Patch <--> Sheet")
                end
                
                imagePatches = [obj.Features(ii).ImagePatch];
                
                imageSheets = cellfun( @(x) strsplit(x, "_WP"), [obj.Features(ii).Sheet], 'UniformOutput', false );
                imageSheets = cellfun( @(x) x(1), imageSheets, 'UniformOutput', false );
                
                indx = ismember([imageSheets{:}], list);
                
                imagePatches = imagePatches(indx);
                
                if size(imagePatches,2) < 2, continue, end
                
                path = fullfile(dir, sprintf('feature%d',count));
                count = count + 1;
                
                if exist(path,'dir') || mkdir(path)
                    
                    iter = 1;
                    for jj = 1: max(num_images)
                        
                        image = imagePatches{iter};
                        image = imresize(image, 0.5, 'method', 'nearest');
                        imwrite(image, fullfile(path, sprintf("instance%d.jpeg",jj)))
                        
                        iter = iter + 1;
                        if iter > size(imagePatches,2), iter = 1; end
                        
                    end
                end
            end
        end
        
    end
    
    
    methods (Access = private)
        % Proigate Features between images
        function obj = Propigate_Features(obj, feature, sheet_)
            
            if nargin > 1
                start = feature;
                stop  = feature;
                sheet = sheet_;
            else
                start = 1;
                stop = size(obj.Features,2);
                sheet = NaN;
            end
            
            
            for numb = start:stop
                
                feature = obj.Features(numb).Boundary;
                
                s= size(obj.Images,2);
                
                count = size(obj.Features(numb).Sheet,1) +1;
                
                for jj = 1:s
                    
                    if jj == sheet, continue, end
                    
                    XX = obj.Images(jj).XX;
                    YY = obj.Images(jj).YY;
                    
                    in = inpolygon(XX,YY, feature(:,1), feature(:,2));
                    in = in  & ~obj.Images(jj).Trough;
                    
                    if any(in(:))
                        
                        obj.Images(jj).Mask(in) = numb;
                        
                        [m,n] = size(obj.Images(jj).Image);
                        
                        mask = zeros(m,n);
                        mask(obj.Images(jj).Mask == numb) = 1;
                        
                        bounds = bwboundaries(mask,'holes');
                        bounds = bounds{1};
                        
                        c = round(mean(bounds(:,1)));
                        r = round(mean(bounds(:,2)));
                        
                        window = SlidingWindow(210, [m,n]);
                        
                        [rows, cols] = window.D2([r,c]);
                        
                        if isempty(rows), continue, end
                        
                        image = obj.Images(jj).Image(rows, cols);
                        
                        %                         figure
                        %                         imshow(image)
                        
                        indx = ismember(obj.Features(numb).Sheet, obj.RawData(jj).log);
                        
                        if any(indx)
                            obj.Features(numb).ImagePatch{indx} = image;
                        else
                            obj.Features(numb).ImagePatch{count} = image;
                            obj.Features(numb).Sheet{count} = obj.RawData(jj).log;
                            count = count +1;
                        end
                        
                    end
                end
            end
            
        end
        
        
        % Reset Image Masks to zeros
        function obj = ResetMasks(obj)
            
            for ii = 1: size(obj.Images,2)
                obj.Images(ii).Mask = zeros(size(obj.Images(ii).Mask));
            end
        end
        
    end
    
    
end


%% Suplimentraty Functions
function xy = RotZ(x, R)

xy = zeros(size(x));

for ii = 1:size(x,1)
    xy(ii,:) = ( R * [x(ii,1);x(ii,2)] )';
    
end
end


% Gepmetric matching Probabilty
function probs = prob_from_cov(mu1, cov1, mu2, cov2)

probs = (mvnpdf(mu1,mu2,cov2)/mvnpdf(mu2,mu2,cov2)) * (mvnpdf(mu2,mu1,cov1)/mvnpdf(mu1,mu2,cov2));

if probs > 0.00001, disp(probs), end

end


% Find overlapping covariance ellipses
function indx = Covariance_mathcing(feature, proposed, plot_io)

indx = false(1, max(size(proposed)));

a1 = ErrorEllipse(feature.xyz(1:2), 1.3 * feature.error, 5);

if plot_io
    fig = figure;
    plot(a1(1,:), a1(2,:), 'k')
    hold on
end

for ii = 1 : max(size(proposed))
    
    a2 = ErrorEllipse(proposed(ii).xyz(1:2), 1.3 * proposed(ii).error, 5);
    in = inpolygon( a1(1,:), a1(2,:), a2(1,:), a2(2,:) );
    
    if any(in), indx(ii) = true; end
    
    if plot_io  && any(in)
        plot(a2(1,:), a2(2,:), 'g')
    elseif plot_io  && ~any(in)
        plot(a2(1,:), a2(2,:), 'r')
    end
    
end

if plot_io
    pause(0.5)
    close(fig)
end

end




