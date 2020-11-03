% Riptide AcommsHandler output parser
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% August 28, 2020



%% Read in pAcommsHandler Logs
function Acomms = Riptide_pAcommsHandler2Mat(file)

inpfmt = 'yyyyMMddHHmmss.SSSSSS';
disfmt = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';

[~, fileName, ext] = fileparts(file);

fprintf(' * Importing Riptide pAcomms Log: %s\n', [fileName,ext]);




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
    Acomms = struct('log','', 'modemID',[], 'msg', []);
    return
    
end

[~,f_name,~] = fileparts(file);

% ____ Get Modem ID ____
temp    = dataArray{3};
place   = strncmp(temp, ': modem_id:',        11);
modemID = temp(place);
modemID = strsplit(modemID{1}, ':');


% ____ Get rid of empty cells _____
dataArray(1) = [];
dataArray(2) = [];

% Times = dataArray{1};
Data  = dataArray{3};


% ____ Parse the remaning data ___
in = strncmp( Data, ': msg: ', 7);

msgs = NaT(sum(in),1, 'Format', disfmt);

count = 1;
for data = Data(in)'
    msg = strsplit(data{1}, {'"',';'});
    msgs(count) = datetime(str2num(msg{2}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', disfmt);
    count = count + 1;
end

dtfns = DateTime_Funs;

msgs = dtfns.Unique(msgs);

% ____ Create Struct _____
Acomms.log         = f_name;
Acomms.modemID     = str2double(modemID{3});
Acomms.msg.msg     = msgs;
Acomms.msg.rx_time = NaT(size(msgs,1), 4, 'Format', disfmt);
Acomms.msg.tx_time = NaT(size(msgs,1), 4, 'Format', disfmt);
Acomms.msg.tof     = NaT(size(msgs,1), 4, 'Format', disfmt);

if isempty(msgs)
    fprintf("\tNo Acomms messages found\n")
    return
end

% ____ Find Interesting info ____
in = ... strncmp( Data, ': D: $CARXP', 11) | ...                                % The modem sends this message when it acoustically detects the start of a packet.
    strncmp( Data, ': D: $CAXST', 11) | ...                                % Communication cycle transmit statistics
    strncmp( Data, ': D: $CACST', 11) | ...                                % Communication cycle receive statistics
    strncmp( Data, ': msg: ',      7)      | ...
    strncmp( Data, ':   source:', 11) | ...
    strncmp( Data, ':   dest:',    9) | ... 
    strncmp( Data, ': D2: Received message: [[ModemTransmission]] src:',  50) | ...
    strncmp( Data, ': D2: Received message: [[ModemTransmission]] time:', 51) | ...
    strcmp(  Data, ': D2: Received trigger: ACOMMS_XMIT');



% ___ Group Communication Packets ____
% disp(Data(in))

data = Data(in);
indx1 = find ( strncmp( data, ': D: $CACST', 11) | ...
              strcmp(  data, ': D2: Received trigger: ACOMMS_XMIT'));

indx2 = [indx1(2:end) - 1; numel(data)];
        
packets{numel(indx1)} = '';
last_packet = {};
ii = 1;
for idx = [indx1';indx2']
    
    packet = data(idx(1):idx(2));
    
    if numel(packet) == 1 && strncmp( packet, ': D: $CACST', 11)
    
        last_packet = packet;
        
    else
        packets{ii} = [last_packet;packet];
        last_packet = {};
        ii = ii + 1;
        
%         disp(packets{ii-1})
    end
    
end

packets(ii:end) = [];




% ____ Extract Data from Packets ____
for packet_ = packets
    
    packet = packet_{1};
    
    idx = find(strncmp( packet, ': msg: ', 7));
    if any(idx)
        msg = mgs2dt( packet{idx(1)}, disfmt);
        [count, ~] = dtfns.Ismember(msgs, msg);
        
    elseif any( strncmp( packet, ': D2: Received message: [[ModemTransmission]] time:', 51) )
        last_packet = packet;
        continue
    else
        disp('No message or time??')
        disp(packet)
        continue
    end
    
    
    % Sending message
    if any( strcmp(  packet, ': D2: Received trigger: ACOMMS_XMIT') )
        
        idx = strncmp( packet, ': D: $CAXST', 11);
        if any(idx)
            dt = caxst2dt(packet(idx), inpfmt, disfmt);
            Acomms.msg.tx_time(count,1:numel(dt)) = dt;
        end
        
        
    % Reciving mesage
    elseif any( strncmp( packet, ': D: $CACST', 11) )
        
        packet = [last_packet; packet];
        last_packet = {};
        
        idx = strncmp( packet, ': D: $CACST', 11);
        if any(idx)
            dt = cacst2dt(packet(idx), inpfmt, disfmt);
            Acomms.msg.rx_time(count,1:numel(dt)) = dt;
        end
        
        
    % Error State
    else
        disp('Random Ass Packet')
        disp(packet)
        
    end
    
    
end


end

%% _____ Old Stuff _____

% count = 1;
% type = '';
% for data = Data(in)'
%    
% %     disp(data)
%     
%    if strncmp( data, ': msg: ', 7)                                          % Message Recived 
%        
%         msg = mgs2dt(data{1}, disfmt);
%         count = find( dtfns.Isequal(msg, msgs) );
% %         disp(count)
%         
%         if strcmp(type, 'rx')                   
%             Acomms.msg.rx_time(count) = rx_time;                            % Rx time comes before the message
%         end
%      
%         
%    elseif strncmp( data, ': D: $CAXST', 11)                                 % Communication cycle transmit statistics
%         type = 'tx';
%         timeStr = strsplit(data{1}, ',');
%         Acomms.msg.tx_time(count) = datetime([timeStr{3},timeStr{4}], 'InputFormat', inpfmt, 'Format', disfmt);
%         
%         
%         
%     elseif strncmp( data, ': D: $CACST', 11)                                % Communication cycle receive statistics
%         type = 'rx';
%         timeStr = strsplit(data{1}, ',');
%         rx_time = datetime(timeStr{4}, 'InputFormat', inpfmt, 'Format', disfmt);  % Transmit message comes befor the time
%         
% 
%  
%         
%     elseif strcmp( data, ': D2: Received trigger: ACOMMS_XMIT')             % About to transmit a message
%         type = 'tx';
%         
%         
%     elseif strncmp( data, ': D2: Received message: [[ModemTransmission]] src:', 50) % Just Recived a message
%         type = 'rx';
%         
%         
%     else                                                                    % Unknow Condition encounterd
%         fprintf('\n\n\n!!! %s !!!\n\n\n', data{1})
%         
%     end
%     
%     
%     
% end
% 
% 
% end






%% Other Functions
function dt = mgs2dt(msg, disfmt)

msg = strsplit(msg, {'"',';'});
dt = datetime(str2num(msg{2}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', disfmt);

end



function dt = caxst2dt(msg, inpfmt, disfmt)

for ii = numel(msg): -1: 1
    timeStr = strsplit(msg{ii}, ',');
    dt(ii) = datetime([timeStr{3},timeStr{4}], 'InputFormat', inpfmt, 'Format', disfmt);
end

end



function dt = cacst2dt(msg, inpfmt, disfmt)

for ii = numel(msg): -1: 1
    timeStr = strsplit(msg{ii}, ',');
    dt(ii) = datetime(timeStr{4}, 'InputFormat', inpfmt, 'Format', disfmt);
end

end

%% Old stuff
% in1 = strncmp( Data, ': D2: Received message: [[ModemTransmission]] src:', 50);
% in2 = strcmp(  Data, ': Successfully decoded message of type: goby.moos.protobuf.RiptideMsg');
%
% if sum(in1) ~= sum(in2)
%     warning( 'Message Indecies do not match!!')
%     return
% end
%
%
% format = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';
%
% fields       = { 'timeStamp', 'timeOfFlight',             'msg',    'src',         'dest' };
% messageTypes = { ':   time:', ':   one_way_travel_time:', ': msg:', ':   source:', ':   dest:' };
% sz           = {  9,           24,                         6,        11,            9 };
%
% count = sum(in1);
% for indx = [find( strncmp( Data, ': D2: Received message: [[ModemTransmission]] src:', 50) )' ;
%             find(  strcmp( Data, ': Successfully decoded message of type: goby.moos.protobuf.RiptideMsg') )']
%
%     data = Data(indx(1):indx(2));
%
%     Acomms.msg(count).timeOfFlight = NaN(1,3);
%
%     for ii = {fields{:}; messageTypes{:}; sz{:}}
%
%         in = strncmp(data, ii{2}, ii{3});
%
%         aa = sum(in);
%         for msg = cellfun( @strsplit, data(in), 'UniformOutput',false)'
%
%             str = msg{1};
%
%             switch ii{1}
%                 case 'msg'
%
%                     b = strsplit(str{3},{';','"'});
%                     d = datetime(str2num(b{2}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', format);
%
%                 case 'timeStamp'
%                     d = datetime(str2num(str{3}), 'ConvertFrom','epochtime','TicksPerSecond',1e6,'Format', format);
%                 otherwise
%                     d = str2double(str{3});
%             end
%
%             Acomms.msg(count).(ii{1})(aa) = d;
%             aa = aa-1;
%
%         end
%
%     end
%     count = count -1;
% end
%
% end





