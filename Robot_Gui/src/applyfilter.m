% Apply filter to to data

function Data = applyfilter(Data, field, filter, filterParameters)

dataInfo = class( Data );

switch dataInfo
    
    case 'struct'
        
        for jj = 1:length(Data)
            
            value = Data(jj).(field);
            
            in = filter(value(:,1),value(:,2),filterParameters(:,2),filterParameters(:,1));
            
            fields = fieldnames(Data);
            
            for ii = 1:length(fields)
                
                if any(strcmp(fields{ii},{'Header','FilteredOut','Log'}))
                    continue
                end
                
                value = Data(jj).(fields{ii});
                
                Data(jj).(fields{ii}) = value(in,:);
                
            end
        end
        
        
    otherwise
        disp('This Object class is not supported:')
        disp(dataInfo)
        
end
end