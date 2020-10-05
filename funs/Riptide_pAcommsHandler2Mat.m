% Riptide AcommsHandler output parser
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% August 28, 2020



function Acomms = Riptide_pAcommsHandler2Mat(file)

[~, fileName, ext] = fileparts(file);

fprintf(' * Importing Riptide pAcomms Log: %s\n', [fileName,ext]);


Acomms = struct('log','', 'modemID',[], 'msg', []);


% ____ Read File ____
try
    formatSpec = '%q%{yyyy-MMM-dd HH:mm:ss.SSSSSS}D%q%q%[^\n\r]';
    fileID = fopen(file,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', {'[',']','{','}'}, 'ReturnOnError', false);
    fclose(fileID);
    
catch exception                 % File Cannot be read --> Skip it and return empty struct
    warning off backtrace
    warning('Error Reading log file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSuspect bad log file, skipping it \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
    warning on backtrace
    return
    
end

[~,f_name,~] = fileparts(file);
Acomms.log = f_name;

% ____ Get Modem ID ____
temp    = dataArray{3};
place   = strncmp(temp, ': modem_id:',        11);
modemID = temp(place);
modemID = strsplit(modemID{1}, ':');

Acomms.modemID = str2double(modemID{3});


% ____ Get rid of unwanted data _____
dataArray(1) = [];
dataArray(2) = [];

Data   = dataArray{3};

% ____ Parse the remaning data ___

in1 = strncmp( Data, ': D2: Received message: [[ModemTransmission]] src:', 50);
in2 = strcmp(  Data, ': Successfully decoded message of type: goby.moos.protobuf.RiptideMsg');


if sum(in1) ~= sum(in2)
    warning( 'Message Indecies do not match!!')
    return 
end


format = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';

fields       = { 'timeStamp', 'timeOfFlight',             'msg',    'src',         'dest' };
messageTypes = { ':   time:', ':   one_way_travel_time:', ': msg:', ':   source:', ':   dest:' };
sz           = {  9,           24,                         6,        11,            9 };

count = sum(in1);
for indx = [find( strncmp( Data, ': D2: Received message: [[ModemTransmission]] src:', 50) )' ;
            find(  strcmp( Data, ': Successfully decoded message of type: goby.moos.protobuf.RiptideMsg') )']
    
    data = Data(indx(1):indx(2));
    
    Acomms.msg(count).timeOfFlight = NaN(1,3);
    
    for ii = {fields{:}; messageTypes{:}; sz{:}}
        
        in = strncmp(data, ii{2}, ii{3});
        
        aa = sum(in);
        for msg = cellfun( @strsplit, data(in), 'UniformOutput',false)'
            
            str = msg{1};
            
            switch ii{1}
                case 'msg'
                    
                    b = strsplit(str{3},{';','"'});
                    d = datetime(str2num(b{2}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', format);
                    
                case 'timeStamp'
                    d = datetime(str2num(str{3}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', format);
                otherwise
                    d = str2double(str{3});
            end
            
            Acomms.msg(count).(ii{1})(aa) = d;
            aa = aa-1;
            
        end
        
    end
    count = count -1;
end

end
