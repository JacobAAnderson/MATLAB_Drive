classdef OceanData
    
    properties
        Data
        LAT
        LON
    end
    
    methods
        
        % Constructor
        function obj = OceanData, warning on backtrace; end
        
        
        % Load Data From a file
        function obj = AddData(obj, param, file)
            
            [~,~,ext] = fileparts(file);
            
            switch ext
                
                % ___ Net CD Files ________________________________________
                case '.nc'
                    
                    info   = ncinfo(file);
                    varnames = {info.Variables(:).Name};
                    
                    ncid = netcdf.open(file);
                    
                    for ii = 1: length(varnames)
                        
                        if strcmp( varnames{ii}, 'palette'), continue, end
                        
                        varid = netcdf.inqVarID(ncid,varnames{ii});
                        data  = netcdf.getVar(ncid,varid);
                        
                        if ~isa(data, 'single') && ~isa(data, 'double')
                            data = single(data);
                        end
                        
                        obj.Data.(param).(varnames{ii}) = data;
                        
                    end
                    
                    if isfield( info, 'Variables') && isfield(info.Variables, 'FillValue')
                        obj.Data.(param).fillvalue = [info.Variables.FillValue]';
                    end
                    
                    netcdf.close(ncid)
                    
                    
                % ___ Hierarchical Data Format files -_____________________
                case '.hdf'
                    S = hdfinfo(file);
                    
                    if isfield(S, 'SDS') && isfield(S.SDS, 'Name')
                        
                        name = S.SDS.Name;
                        
                        obj.Data.(param).(name) = hdfread(file, name);
                    end
                    
                    
                    if isfield(S, 'Attributes') && isfield(S.Attributes, 'Name')
                        names = {S.Attributes.Name};
                        starIdx = strcmpi(names, 'Start Time String');
                        stopIdx = strcmpi(names, 'Stop Time String');
                        
                        infmt = 'MM/dd/yyyy HH:mm:ss';
                        
                        obj.Data.(param).startime = datetime(S.Attributes(starIdx).Value, 'InputFormat',infmt);
                        obj.Data.(param).stoptime = datetime(S.Attributes(stopIdx).Value, 'InputFormat',infmt);
                    else
                        warning('No start or stop attributes present for %s', param)
                    end
                    
                    if isfield(S, 'SDS') && isfield(S.SDS, 'Attributes')
                        obj.Data.(param).fillvalue = S.SDS.Attributes(16).Value;
                        obj.Data.(param).limets = S.SDS.Attributes(7).Value;
                    else
                        warning('No fill values or limets present for %s', param)
                    end
                    
                    
                % ___ Matlab Data File ____________________________________    
                case '.mat'
                    var = load(file);
                    
                    fds = fields(var);
                    if numel(fds) == 1
                        obj.Data.(param).(param) = var.(fds{1});
                    else
                        for ii = 1: numel(fds)
                            obj.Data.(param).(fds{ii}) = var.(fds{ii});
                        end
                    end
                    
                    obj.Data.(param).fillvalue = Inf;                       % Assume a fill value of Inf
                    
                % ___ Wrong file type _____________________________________    
                otherwise
                    warning('Unsupported file Type')
            end
            
        end
        
        
        % Apply Function to Modify Data
        function obj = Modify(obj, param, func)
            
            if ~isfield(obj.Data, param)
                warning("Unreckognized Parameter %s", param)
                return
            end
            
            obj.Data.(param).(param) = func(obj.Data.(param).(param));
            
        end
        
        
        % Construct World Coordinate system
        function obj = World_Coordinates(obj, params)
            
            n = numel(params);                                              % Number of Parameters to be used
            
            lat{n} = zeros(10);                                             % Peralocate Some memory
            lon{n} = zeros(10);
            
            tf_lon = true;
            tf_lat = true;
            for ii = 1: n                                                   % Collect the avalable cordinate systems
                if isfield(obj.Data, params{ii}) && isfield(obj.Data.(params{ii}),'lat') && isfield(obj.Data.(params{ii}),'lon')
                    
                    lat{ii} = obj.Data.(params{ii}).lat;                    % Collect varavles incase further manipulation is needed
                    lon{ii} = obj.Data.(params{ii}).lon;
                    
                    tf_lat = isequal(lat{1}, lat{ii}) & tf_lat;             % Check to see if all the arrays are the same
                    tf_lon = isequal(lon{1}, lon{ii}) & tf_lon;
                    
                else
                    fprintf('Field %s Does Not Contain Lat Lon Data', params{ii})
                end
                
                
            end
            
            
            if tf_lon && tf_lat &&  ismatrix(lat{1}) && min(size(lat{1})) == 1 && min(size(lon{1})) == 1 % All the reference systems match
                [obj.LAT, obj.LON] = meshgrid(double(lat{1}), double(lon{1})); 
            else  
            	warning("I don't know wat to do!!!!")
            end
        end
        
        
        % Fit Data to world coordinate system
        function obj = Fit2World(obj, params)
            
            if isempty(obj.LON)
                warning('World Coordinate system has not been established')
                return
            end
            
            
            n = numel(params);
            
            newData(:,:,12) = obj.LAT;
            
            for ii = 1:n
                
                data = obj.Data.(params{ii}).(params{ii});
                fill = obj.Data.(params{ii}).fillvalue(1);
                
                if any( size(data) ~= size(obj.LON))
                    
                    switch ndims(data)
                        case 2, S = 1;
                        case 3, S = size(data,3);
                        otherwise, warning('The dimentions of %s are not suitable for conversion')
                                   continue
                    end
                    
                    [r,c] = size(data);
                    
                    x = linspace(-180, 180, c); %(1:c) * 360/c -180;
                    y = linspace(-90, 90, r);  %(1:r) * 180/c -90;
                    
                    [X,Y] = meshgrid(x,y);
                    
                    newData(:,:,S) = obj.LAT;
                    
                    for s = 1:S
                        thisData = double(data(:,:,s));
                        thisData(thisData == fill) = NaN;
                        newData(:,:,s) = interp2(X,Y, thisData, obj.LON, obj.LAT, 'natural');
                    end
                    
                end
                
                obj.Data.(params{ii}).(params{ii}) = newData(:,:,1:s);
            end
            
            
            
        end
        
        
        % Calculate daylength
        function obj = Daylength(obj, param)
            lat = obj.LAT(1,:);
            yd = day(obj.Data.(param).startime,'dayofyear');
            
            obj.Data.dl.dl = daylength(yd,lat);
           
        end
        
        
        % Access Data Values
        function data = Value(obj, param, scale)
            
            if ~isfield(obj.Data, param)
                warning("Unreckognized Parameter %s", param)
                return
            end
            
            data = obj.Data.(param).(param);
            
            if nargin > 2, data = imresize(data,scale); end
            
        end
        
        
        % Plot Data
        function fig = Plot(obj, param, fig)
            
            if ~isfield(obj.Data, param)
                warning("Unreckognized Parameter %s", param)
                return
            end
            
            if nargin > 2, figure(fig);
            else, fig = figure('name', param, 'numbertitle', 'off');
            end
            
            J = obj.Data.(param).(param);
            fill = obj.Data.(param).fillvalue(1);
            
            J(J == fill) = NaN;
            
            if ~isempty(obj.LON)
                
                scale = 0.25;
                LON_ = imresize(obj.LON,scale);
                LAT_ = imresize(obj.LAT,scale);
                J   = imresize(J, scale);
                
                mesh(LON_, LAT_, J)
                
                
            elseif isfield(obj.Data.(param),'lat') && isfield(obj.Data.(param),'lon')
                
                [LAT_, LON_] = meshgrid(obj.Data.(param).lat, obj.Data.(param).lon);
                
                scale = 0.25;
                LON_ = imresize(LON_,scale);
                LAT_ = imresize(LAT_,scale);
                J   = imresize(J, scale);
                
                mesh(LON_, LAT_, J)
                
            else
                mesh(J)
            end
            
            title(param, 'Interpreter','none')
            colorbar
            view(0,90)
            
        end
        
    end
end