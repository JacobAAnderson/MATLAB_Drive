% Model Acoustic Communication Latency
% Jacob Anderson
% Sept 22, 2020


function [mu, sig] = Model_AcommsLatency( RT )

dtfs = DateTime_Funs;                                                       % Date Time Functions
disfmt = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';

fprintf('\n\nLatency in acoustic Communications:\n\n')

% ---- Colect Acoms data --------------------------------------------------
v1 = cat(1,RT(1).rawAcomms.vehicle);                                        % Vehicle location
v2 = cat(1,RT(2).rawAcomms.vehicle);

sz(1) = size(v1,1);                                                         % Array Sizes
sz(2) = size(v2,1);
sz(3) = max(sz(1), sz(2));

lat = NaN(sz(3),2);                                                         % Vehicle Latitudes
lat(1:sz(1),1) = v1(:,1);
lat(1:sz(2),2) = v2(:,1);

lon = NaN(sz(3),2);                                                         % Vehicle Longitudes
lon(1:sz(1),1) = v1(:,2);
lon(1:sz(2),2) = v2(:,2);

tx  = NaT(sz(3),2, 'Format', disfmt);
tx(1:sz(1),1) = cat(1,RT(1).rawAcomms.sent_msg);                            % Sent Messages
tx(1:sz(2),2) = cat(1,RT(2).rawAcomms.sent_msg);

rx  = NaT(sz(3),2, 'Format', disfmt);
rx(1:sz(1),1) = cat(1,RT(1).rawAcomms.recived_msg);                         % Sent Messages
rx(1:sz(2),2) = cat(1,RT(2).rawAcomms.recived_msg);

ts = NaT(sz(3),2, 'Format', disfmt);
ts(1:sz(1),1) = cat(1,RT(1).rawAcomms.timeStamp);                           % Time Stamps
ts(1:sz(2),2) = cat(1,RT(2).rawAcomms.timeStamp);

tof = NaT(sz(3),2, 'Format', disfmt);
tof(1:sz(1),1) = cat(1,RT(1).rawAcomms.tof);                                % Time of flight
tof(1:sz(2),2) = cat(1,RT(2).rawAcomms.tof);

a = numel(RT);

mu(a)  = seconds;
sig(a) = seconds;

for a = [1:a; a:-1:1]
    
    s = a(1);                                                               % Sending  Vehicle
    r = a(2);                                                               % Reciving Vehicle

    % ---- Corelate Transmited and recived messsages ----------------------
    [tf, ind] = dtfs.Ismember( tx(:,s), rx(:,r) );                          % Vehicle A recived communications from vehicle B
     
    ind(ind == 0) = [];
    
    set = [tx(tf,s), rx(ind,r),];                                           % Messages from vehicle A recived by vehicle B
    
    if ~dtfs.Isequal(set(:,1), set(:,2)), warning("Set Indexing is wrong"), end
    
    
    % ---- Calculate One Way Time of Flight from time stamp ---------------
%     t_sent   = tx(tf,s);                                                    % Time message was sent, corrected for clock error
%     t_recive = ts(ind,r);                                                   % Time message was recived
    
    t_sent   = tof(tf,s);                                                   % Time message was sent, corrected for clock error
    t_recive = tof(ind,r);                                                  % Time message was recived
    
    owtof = t_recive - t_sent;                                              % One way time of flight
%     owtof(isnan(owtof)) = [];                                               % Get rid of NaTs
%     owtof = owtof;
        
    
    % ---- Calulate distance between vehicles and true time of flight -----
    d = vdist(lat(tf,s),lon(tf,s),lat(ind,r),lon(ind,r));                   % Distance between vehicles
    
    t = seconds(d/1475);                                                    % Time of flight based on speed of sound in water --> 1475 [m/s]
    
    
%     %% !!! Compare pAcomms TOF HERE !!!
%     pTOF = cat(1, RT(r).rawAcomms.tof);
%     pTOF = pTOF(ind,:);
    
    
    % ---- Calculate the acomms latency -----
    err = seconds(t - owtof);                                                          % Latency between sending the message and reciving the message
    
    mu(s)  = seconds( nanmean(err));                                                     % Average latency
    sig(s) = seconds( nanstd(err));                                                      % Standard deviation in latency

    fprintf('Vehicle %s to %s\n', RT(s).name, RT(r).name)
    fprintf('Ave: %s,\tSTD: %s\n',mu(s), sig(s))
    fprintf('Population: %d\n\n', size(set,1))
    
end

fprintf('\n\n')

end




