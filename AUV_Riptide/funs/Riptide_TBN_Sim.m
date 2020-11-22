% Riptide DEC-TBN Simulation
% Jacob Anderson
% RDML, OSU, Corvallis, OR.
% November 6, 2020




function [tbnPaths, sim_times, RT] = Riptide_TBN_Sim(RT, timeLine, sos, mu, dec_tbn)

dtfs = DateTime_Funs;                                                       % Date Time Functions

count = [0,0];                                                              % Count number of Acoustic communications

comms_count = [0,0];

sim_times = NaT( size(timeLine,1), 2);


% ---- Run through data logs doing TBN -----

fprintf('\n==== Performing %s =========================\n', dec_tbn)
    
for time = timeLine'                                                        % Run the timeline
   
    for v = [1:numel(RT); numel(RT):1]                                      % Switch between vehicles
             
        for mem = {'data', 'acomms'}                                                  % Switch between data sources
    
            for in = find( dtfs.Ismember( RT(v).(mem{1}).gpsDate, time) )'  % Find which data it is time for
                
                % ---- Process Navigation data ----
                if strcmp(mem{1},'data')
                    
                    % --- Particle Filter Update ---
                    dt      =  RT(v).data.timeStep(in);                     % Get time step
                    speed   =  RT(v).filteredData.speed(in);                % Use the filtered speed, heading and altitude 
                    heading =  RT(v).filteredData.heading(in);
                    
                    alt     = -RT(v).filteredData.altitude(in);             % Use Filtered Data
                    offset  =  RT(v).filteredData.alt_offset(in,:);
                    
%                     alt     = -RT(v).orical.alt(in);                        % Feed in perfect Altimiter readings
                     
                    if alt == 0 || isnan(alt)                               % Check for valid altimiter readings
                        RT(v).tbn = RT(v).tbn.Update(speed, dt, heading);
                    else
                        RT(v).tbn = RT(v).tbn.Update(speed, dt, heading, alt, offset);
                    end
                   
                    count(v) = count(v) + 1;
                    
                    sim_times( count(v), v ) = time;
                    
                end
                
                
                
                % ---- Process Recived Acomms messages ----
                if strcmp(mem{1},'acomms') && ~isnat(RT(v).acomms.recived_msg(in)) && RT(v).acomms.policy.(dec_tbn)(in)
                    
                    if v == 1, s=2;                                         % Get index of the vehicle that sent the message
                    else, s = 1;
                    end
                    
                    owtof = RT(v).acomms.tof(in) - mu(v);                   % One way time of flight
                    dist  = seconds(owtof) * sos;                           % Compute Distance
                    x     = RT(s).tbn.X(1:2)';
                    sig   = RT(s).tbn.cov(1:2,1:2);
                    
                    fprintf("[%s]  TBN - Acomms, owtof: %1.8s, dist: %5.2f [m]\n",RT(v).name, owtof, dist) 
                    
                    RT(v).tbn = RT(v).tbn.Acoms(dist, x, sig);
                    
                    comms_count(v) = comms_count(v) + 1;
                    
                    
                end
                
                
            end
        end
    end
end


fprintf('\nNumber of Communications Recived By\n')
fprintf('%s, %d\n%s, %d\n\n', RT(1).name, comms_count(1), RT(2).name, comms_count(2) )

% --- Copy over TBN path to RT and conver to lat_lon ---
    
sz = max( size(RT(1).tbn.path,1), size(RT(2).tbn.path,1) ); 

tbnPaths = zeros(sz,2,numel(RT));

for v = 1: numel(RT)                                                
    sz = size(RT(v).tbn.path,1);
    tbnPaths(1:sz,:,v) = RT(v).tbn.path(:,1:2);
    
end


end




