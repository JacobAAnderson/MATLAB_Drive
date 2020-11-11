% Riptide Acomms
% Jacob Anderson
% DRML, OSU, Corvallis, OR.
% Sept 22, 2020



classdef Riptide_Acomms
    
    
    methods
        
        % Model Acoustic Communication Latency
        function [mu, sig, OWTOF] = Model_AcommsLatency(~, RT )
            
            dtfs = DateTime_Funs;                                           % Date Time Functions
            disfmt = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';
            
            fprintf('\n\nLatency in acoustic Communications:\n\n')
            
            % ---- Colect Acoms data --------------------------------------------------
            v1 = cat(1,RT(1).rawAcomms.vehicle);                            % Vehicle location
            v2 = cat(1,RT(2).rawAcomms.vehicle);
            
            sz(1) = size(v1,1);                                             % Array Sizes
            sz(2) = size(v2,1);
            sz(3) = max(sz(1), sz(2));
            
            lat = NaN(sz(3),2);                                             % Vehicle Latitudes
            lat(1:sz(1),1) = v1(:,1);
            lat(1:sz(2),2) = v2(:,1);
            
            lon = NaN(sz(3),2);                                             % Vehicle Longitudes
            lon(1:sz(1),1) = v1(:,2);
            lon(1:sz(2),2) = v2(:,2);
            
            tx  = NaT(sz(3),2, 'Format', disfmt);
            tx(1:sz(1),1) = cat(1,RT(1).rawAcomms.sent_msg);                % Sent Messages
            tx(1:sz(2),2) = cat(1,RT(2).rawAcomms.sent_msg);
            
            rx  = NaT(sz(3),2, 'Format', disfmt);
            rx(1:sz(1),1) = cat(1,RT(1).rawAcomms.recived_msg);             % Sent Messages
            rx(1:sz(2),2) = cat(1,RT(2).rawAcomms.recived_msg);
            
            % ts = NaT(sz(3),2, 'Format', disfmt);
            % ts(1:sz(1),1) = cat(1,RT(1).rawAcomms.timeStamp);             % Time Stamps
            % ts(1:sz(2),2) = cat(1,RT(2).rawAcomms.timeStamp);
            
            tof(1).t = NaT(sz(3),4, 'Format', disfmt);
            tof(2).t = NaT(sz(3),4, 'Format', disfmt);
            tof(1).t(1:sz(1),:) = cat(1,RT(1).rawAcomms.pAcomms_tof);       % Time of flight
            tof(2).t(1:sz(2),:) = cat(1,RT(2).rawAcomms.pAcomms_tof);
            
            a = numel(RT);
            
            mu(a)  = seconds;
            sig(a) = seconds;
            
            if nargout ==3
                OWTOF(a).time = seconds;
            end
            
            for a = [1:a; a:-1:1]
                
                s = a(1);                                                   % Sending  Vehicle
                r = a(2);                                                   % Reciving Vehicle
                
                % ---- Corelate Transmited and recived messsages ----------------------
                [tf, ind] = dtfs.Ismember( tx(:,s), rx(:,r) );              % Vehicle A recived communications from vehicle B
                
                ind(ind == 0) = [];
                
                set = [tx(tf,s), rx(ind,r),];                               % Messages from vehicle A recived by vehicle B
                
                dtfs.SetError(set(:,1), set(:,2));
                
                
                % ---- Calculate One Way Time of Flight from time stamp ---------------
                %     t_sent   = tx(tf,s);                                  % Time message was sent, corrected for clock error
                %     t_recive = ts(ind,r);                                 % Time message was recived
                
                t_sent   = tof(s).t(tf,:);                                  % Time message was sent, corrected for clock error
                t_recive = tof(r).t(ind,:);                                 % Time message was recived
                
                if ~ isequal( isnat(t_sent), isnat(t_recive) )              % Correct for mis matches in data entry
                    disp('Miss match in TOFs')
                    
                    % Shift rows with 3 entries forward by 1 coulmn
                    idx = sum(isnat(t_recive), 2) == 1;
                    
                    t_recive(idx,1) = t_recive(idx,2);
                    t_recive(idx,2) = t_recive(idx,3);
                    
                    
                    % Shift rows with 4 entries forward by 2 coulmn
                    idx = sum(isnat(t_recive), 2) == 0;
                    
                    t_recive(idx,1) = t_recive(idx,3);
                    t_recive(idx,2) = t_recive(idx,4);
                    
                end
                
                
                owtof_ = t_recive - t_sent;                                 % One way time of flight
                
                % ---- Calulate distance between vehicles and true time of flight -----
                d = vdist(lat(tf,s),lon(tf,s),lat(ind,r),lon(ind,r));       % Distance between vehicles
                
                t = seconds(d/1475);                                        % Time of flight based on speed of sound in water --> 1475 [m/s]
                
                
                % ---- Calculate the acomms latency -----
                err = seconds(t - owtof_);                                  % Latency between sending the message and reciving the message
                
                mu(s)  = seconds( nanmean(err(:)) );                        % Average latency
                sig(s) = seconds( nanstd(err(:))  );                        % Standard deviation in latency
                
                fprintf('Vehicle %s to %s\n', RT(s).name, RT(r).name)
                fprintf('Ave: %s,\tSTD: %s\n',mu(s), sig(s))
                fprintf('Population: %d\n\n', size(set,1))
                
                if nargout ==3
                    OWTOF(s).time = nanmean(owtof_, 2);
                    OWTOF(s).msgs = tx(tf,s);
                    OWTOF(s).mu   = mu(s);
                    OWTOF(s).sig  = sig(s);
                    OWTOF(s).gt   = d;
                    OWTOF(s).time.Format = 's';
                end
            end
            
            fprintf('\n\n')
            
        end
        
        
        
        function [mu, sig, d] = Mission_Latency(~, RT )
            
            dtfs = DateTime_Funs;                                           % Date Time Functions
            disfmt = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';
            
            fprintf('\n\nLatency in acoustic Communications:\n\n')
            
            % ---- Colect Acoms data --------------------------------------------------
            v1 = cat(1,RT(1).acomms.vehicle);                               % Vehicle location
            v2 = cat(1,RT(2).acomms.vehicle);
            
            sz(1) = size(v1,1);                                             % Array Sizes
            sz(2) = size(v2,1);
            sz(3) = max(sz(1), sz(2));
            
            lat = NaN(sz(3),2);                                             % Vehicle Latitudes
            lat(1:sz(1),1) = v1(:,1);
            lat(1:sz(2),2) = v2(:,1);
            
            lon = NaN(sz(3),2);                                             % Vehicle Longitudes
            lon(1:sz(1),1) = v1(:,2);
            lon(1:sz(2),2) = v2(:,2);
            
            tx  = NaT(sz(3),2, 'Format', disfmt);
            tx(1:sz(1),1) = cat(1,RT(1).acomms.sent_msg);                   % Sent Messages
            tx(1:sz(2),2) = cat(1,RT(2).acomms.sent_msg);
            
            rx  = NaT(sz(3),2, 'Format', disfmt);
            rx(1:sz(1),1) = cat(1,RT(1).acomms.recived_msg);                % Sent Messages
            rx(1:sz(2),2) = cat(1,RT(2).acomms.recived_msg);
            
            tof(1).t = NaT(sz(3),4, 'Format', disfmt);
            tof(2).t = NaT(sz(3),4, 'Format', disfmt);
            tof(1).t(1:sz(1),:) = cat(1,RT(1).acomms.pAcomms_tof);          % Time of flight
            tof(2).t(1:sz(2),:) = cat(1,RT(2).acomms.pAcomms_tof);
                       
            a = numel(RT);
            
            mu(a)  = seconds;
            sig(a) = seconds;
            
            for a = [1:a; a:-1:1]
                
                s = a(1);                                                   % Sending  Vehicle
                r = a(2);                                                   % Reciving Vehicle
                
                % ---- Corelate Transmited and recived messsages ----------------------
                [tf, ind] = dtfs.Ismember( tx(:,s), rx(:,r) );              % Vehicle A recived communications from vehicle B
                
                ind(ind == 0) = [];
                
                set = [tx(tf,s), rx(ind,r),];                               % Messages from vehicle A recived by vehicle B
                
                dtfs.SetError(set(:,1), set(:,2));
                
                
                % ---- Calculate One Way Time of Flight from time stamp ---------------
                %     t_sent   = tx(tf,s);                                  % Time message was sent, corrected for clock error
                %     t_recive = ts(ind,r);                                 % Time message was recived
                
                t_sent   = tof(s).t(tf,:);                                  % Time message was sent, corrected for clock error
                t_recive = tof(r).t(ind,:);                                 % Time message was recived
                
                if ~ isequal( isnat(t_sent), isnat(t_recive) )              % Correct for mis matches in data entry
                    disp('Miss match in TOFs')
                    
                    % Shift rows with 3 entries forward by 1 coulmn
                    idx = sum(isnat(t_recive), 2) == 1;
                    
                    t_recive(idx,1) = t_recive(idx,2);
                    t_recive(idx,2) = t_recive(idx,3);
                    
                    
                    % Shift rows with 4 entries forward by 2 coulmn
                    idx = sum(isnat(t_recive), 2) == 0;
                    
                    t_recive(idx,1) = t_recive(idx,3);
                    t_recive(idx,2) = t_recive(idx,4);
                    
                end
                
                
                owtof_ = t_recive - t_sent;                                 % One way time of flight
                
                % ---- Calulate distance between vehicles and true time of flight -----
                d = vdist(lat(tf,s),lon(tf,s),lat(ind,r),lon(ind,r));       % Distance between vehicles
                
                t = seconds(d/1475);                                        % Time of flight based on speed of sound in water --> 1475 [m/s]
                
                
                % ---- Calculate the acomms latency -----
                err = seconds(t - owtof_);                                  % Latency between sending the message and reciving the message
                
                mu(s)  = seconds( nanmean(err(:)) );                        % Average latency
                sig(s) = seconds( nanstd(err(:))  );                        % Standard deviation in latency
                
                fprintf('Vehicle %s to %s\n', RT(s).name, RT(r).name)
                fprintf('Ave: %s,\tSTD: %s\n',mu(s), sig(s))
                fprintf('Population: %d\n\n', size(set,1))
                
  
            end
            
            fprintf('\n\n')
            
        end
        
        
        
        function figOut = Plot_AcousticCommunications(~, RT, geotiff)
            
            dtfs = DateTime_Funs;
            
            % --- Extract Data ----
            name1 = RT(1).name;                                                         % Vehicle Name
            path1 = RT(1).data.vehicle(:,1:2);                                          % Complete GPS Track [ lat, lon]
            lat1  = RT(1).acomms.vehicle(:,1);                                          % Location of acomms
            lon1  = RT(1).acomms.vehicle(:,2);
            tx1   = RT(1).acomms.sent_msg;                                              % Sent Message
            rx1   = RT(1).acomms.recived_msg;                                           % Recived Message
            
            
            name2 = RT(2).name;
            path2 = RT(2).data.vehicle(:,1:2);
            lat2  = RT(2).acomms.vehicle(:,1);
            lon2  = RT(2).acomms.vehicle(:,2);
            tx2   = RT(2).acomms.sent_msg;
            rx2   = RT(2).acomms.recived_msg;
            
            
            % --- Corrolate the communications ---
            [Lia,Locb] = dtfs.Ismember(tx1,rx2);
            Locb(Locb == 0) = [];
            
            [Lib,Loca] = dtfs.Ismember(tx2,rx1);
            Loca(Loca == 0) = [];
            
            
            % --- Make Sure the indexing is right ---
            set = [tx1(Lia), rx2(Locb),];
            
            dtfs.SetError(set(:,1), set(:,2))
            
            
            % --- Plot ---
            fig = figure('Name',"Acomms Map", 'numbertitle','off');
            
            if nargin == 3, geotiff.Show; end
            
            hold on
            
            % Plot vehicle GPS tracks
            p1 = plot(path1(:,2), path1(:,1), 'y');
            p2 = plot(path2(:,2), path2(:,1), 'b');
            
            % Plot Acoustic Communications sites
            plot(lon1, lat1, '*y')
            plot(lon2, lat2, '*b')
            
            
            % Plot Communications from vehicle A to B
            for arr = [lon1(Lia)'; lat1(Lia)'; zeros(size(lon1(Lia)))'; lon2(Locb)'; lat2(Locb)'; zeros(size(lon2(Locb)))']
                
                farrow( arr(1), arr(2), arr(3), arr(4), arr(5), arr(6), [0.6314    0.9882    0.3765], 1 );
                
            end
            
                        
            % Plot Communications from vehicle B to A
            for arr = [lon2(Lib)'; lat2(Lib)'; zeros(size(lon2(Lib)))'; lon1(Loca)'; lat1(Loca)'; zeros(size(lon1(Loca)))']
                
                farrow( arr(1), arr(2), arr(3), arr(4), arr(5), arr(6), 'c', 1 );
                
            end
            
            hold off
            
            legend([p1,p2], {name1, name2})
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
            
        end
        
        
    end
    
end


