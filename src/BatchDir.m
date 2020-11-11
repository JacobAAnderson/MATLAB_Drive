

function Data = BatchDir(filePath,fileType,fun)

fprintf('\nBegin Importing Data From: %s\n', filePath)

files = fullfile(filePath,['*.',fileType]);
logfiles = dir(fullfile(files));

fprintf('\n %d  %s Files Found\n',size(logfiles,1),fileType);

if isempty(logfiles)
    Data = [];
    return
end


filePath = cellstr(repmat(filePath, size(logfiles,1),1));                   % Collect file names
files = fullfile(filePath, {logfiles.name}');


% emptyFiles = cat(1, logfiles.bytes) < 950;                                  % Get rid of empty files
% files(emptyFiles) = [];
% 
% if sum(emptyFiles) > 0
%     fprintf('\n ---> %d Empty Files Filtered Out\n',sum(emptyFiles))
% end


A = cellfun( fun , files, 'UniformOutput', false );                         % Process Files


switch class(A{1})
    
    case 'struct'
        
        Data = cat(1, A{:});                                                % Extract data from the cell array
        
        % Sort Data into chronological order ------------------------------------
        if any(isfield(Data,'timeStamp')) 
            
            sortArray(size(Data,1)) = datetime;
            
            for ii = 1: size(Data,1)
                sortArray(ii) = Data(ii).timeStamp(1);
            end
            
            % Get rid of empty bad data entries
            in = ~isnat(sortArray);
            Data = Data(in);
            sortArray = sortArray(in);
            
            [~,I] = sort(sortArray);
            
            Data = Data(I);
        end
        
    otherwise
        
        Data = A;
        
end

fprintf('\nFinished Importing Data \n\n')
end

