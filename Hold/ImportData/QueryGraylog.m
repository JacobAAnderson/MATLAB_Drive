clc
close all
% clear all

% https://logger.fortlewis.edu:12900/api/search/universal/keyword?query=*&keyword=last%2090%20days&fields=ph1
% https://logger.fortlewis.edu:12900/api/search/universal/keyword?query=Data&keyword=last%2090%20days

% url = 'https://logger.fortlewis.edu:12900/api/search/universal/keyword?query=Data&keyword=last%2090%20days&messages=200';


searchUrl = 'https://logger.fortlewis.edu:12900/api/search/universal/keyword?';
query = 'query=Data&keyword=last%2090%20days&messages=200';


options = weboptions('Username','iv0fiugrs9gju83i84547pl9iba1bv8pidue33spktosrjl56v4',...
    'Password','token', ...
    'MediaType','application/json');


data = webread([searchUrl,query], options);
% fprintf('%s',data.built_query)

%% Extract data from messages
clc
[m,~] = size(data.messages);        % Get number of messages

Timestamp(m) = datetime;
Timestep = 10;
GPS     = zeros(m,2);                % Preallocate memory
Battery = zeros(m,1);
Array1  = zeros(m,11);
Array2  = zeros(m,11);

ii = 1;

for jj = 1:m                        % Iterae throught messages and assign values
    
    datum = data.messages(jj).message.full_message;
    parts = strsplit(datum,';');
    
    gps = strsplit( parts{1},',');
    
    if ~isempty(gps{1})           % fill in missing Date info
        day = gps{1};
    end
    
    try
        Timestamp(ii) = datetime([day,' ',gps{2}],'InputFormat', 'ddMMyy HHmmss.SSS');
    catch
        Timestamp(ii) = Timestamp(ii-1) + minutes(10);
    end
    
    GPS(ii,:)     = [str2double(gps(3)), str2double(gps(4))];
    Battery(ii)   = str2double(parts{2});
    Array1(ii,:)  = str2double(strsplit( parts{3}, ','));
    Array2(ii,:)  = str2double(strsplit( parts{4}, ','));
    
    Array1(ii,10) = -Array1(ii,10);
    Array2(ii,10) = -Array2(ii,10);
    
    ii = ii+1;
    
end

Timestamp(ii:m) = [];
GPS(ii:m,:)     = [];
Battery(ii:m)   = [];
Array1(ii:m,:)  = [];
Array2(ii:m,:)  = [];

clear datum;
clear parts
clear gps;
clear day;
clear ii;

Timestamp = DateTimeFilter(Timestamp);
Array1 = LowPassFilter(Array1, 2);
Array1 = LowPassFilter(Array1, 2);

%% Plotting
close all

figure('Name', 'Battery Volage','NumberTitle','off');
plot( Timestamp', Battery)
ylabel('Volts')


figure('Name', 'Node 1 Depth and Temperature','NumberTitle','off');
hold on
subplot(2,1,1)
plot( Timestamp', Array1(:,11));
ylabel('Degrees C')

subplot(2,1,2)
plot( Timestamp', Array1(:,10))
ylabel('Meter')
hold off


figure('Name', 'Node 2 Depth and Temperature','NumberTitle','off');
hold on
subplot(2,1,1)

plot( Timestamp', Array2(:,11))
ylabel('Degrees C')

subplot(2,1,2)
plot( Timestamp', Array2(:,10))
ylabel('Meters')


