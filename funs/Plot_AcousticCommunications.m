% Plot Acoustic Communication
% Jacob Anderson
% Sept 22, 2020


function fig = Plot_AcousticCommunications(RT, geotiff)

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
[Lia,Locb] = ismember(tx1,rx2);                                             
Locb(Locb == 0) = [];

[Lib,Loca] = ismember(tx2,rx1);
Loca(Loca == 0) = [];


% --- Make Sure the indexing is right ---
set = [tx1(Lia), rx2(Locb),];                                           

if ~isequal(set(:,1), set(:,2)), warning("Set Indexing is wrong"), end


% --- Plot ---
fig = figure('Name',"Acomms Map", 'numbertitle','off');

if nargin == 2, geotiff.Show; end

hold on

p1 = plot(path1(:,2), path1(:,1), 'y');
p2 = plot(path2(:,2), path2(:,1), 'b');
plot(lon1, lat1, '*y')
plot(lon2, lat2, '*b')


% Plot Communications from vehicle A to B
for arr = [lon1(Lia)'; lat1(Lia)'; zeros(size(lon1(Lia)))'; lon2(Locb)'; lat2(Locb)'; zeros(size(lon2(Locb)))']

    farrow( arr(1), arr(2), arr(3), arr(4), arr(5), arr(6), [0.6314    0.9882    0.3765], 1 );

end



% Plot Communications from vehicle B to A
d_lat = lat1(Loca) - lat2(Lib);
d_lon = lon1(Loca) - lon2(Lib);
quiver( lon2(Lib), lat2(Lib), d_lon, d_lat, 1)

for arr = [lon2(Lib)'; lat2(Lib)'; zeros(size(lon2(Lib)))'; lon1(Loca)'; lat1(Loca)'; zeros(size(lon1(Loca)))']

    farrow( arr(1), arr(2), arr(3), arr(4), arr(5), arr(6), 'c', 1 );

end

hold off

legend([p1,p2], {name1, name2})

end



