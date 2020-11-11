% Jacob Anderson
% Google Earth KML to matalb data srtucture
% Sept 16, 2020


function Data = GoogleEarthKML2mat(file)

Data = struct('file','','name', '', 'coordinates', []);

[~, fileName, ext] = fileparts(file);

fprintf(' * Reading Google Eath File: %s\n', [fileName,ext]);
            
Data.file = [fileName,ext];    % Add name of the log file to the data structure


% Check that the file path is valid
if ~exist(file,'file')
    fprintf('\tThe File %s Does Not Exist\n', fileName);
    return
end

% Check that the files has data. i.e.. more than just a header
logfile = dir(file);
if (logfile.bytes < 1000 )
    fprintf('\t--> File is empty\n')
    return
end


% Read txt file
try 
    fileID = fopen(file,'r');
    formatSpec = '%[^\n\r]';                                                      % Determin data format based on variable
    dataArray = textscan(fileID, formatSpec,  'ReturnOnError', false); % read the rest of the file
    fclose(fileID);
    dataArray = dataArray{1};
    
catch exception                 % File Cannot be read --> Skip it and return empty struct
    warning off backtrace
    warning('Error Reading log file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSuspect bad log file, skipping it \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
    warning on backtrace
    return
    
end


types = {'<Polygon>','<LineString>'};

count = 1;

Data(100).name = '';
Data(100).coordinates = [];

for range = [find(strcmp( '<Placemark>', dataArray))';
             find(strcmp( '</Placemark>',dataArray))']
            
    
    data = dataArray(range(1):range(2));
    
    tf = ismember(types, data);
    
    if ~any(tf), continue, end                                     % Skip points and unknow info, we only want polygons and paths
    
    
    ind = strncmp('<name>', data, 6);
    
    name = strsplit(data{ind}, {'<','>'});
    
    Data(count).name = name{3};
    
    ind = find(strncmp('<coordinates>', data, 13));
        
    numbs = str2num(data{ind+1});
    
    Data(count).coordinates = reshape(numbs, 3, [])';
        
    fprintf('\tAdding %s: %s\n', name{3}, types{tf})
        
    count = count + 1;
    
end

Data(count:end) = [];                                                       % Get Rid of extra space

end