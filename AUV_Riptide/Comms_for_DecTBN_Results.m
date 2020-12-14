%% Plot Comms planning for dec-TBN results
% Jake Anderson
% Dec 7, 2020



%% Field Trials
clear all
close all
clc
mission1 = [135.0,  87.0, 78.0, 55.0];
mission2 = [ 52.5,  55.7, 66.9, 49.1];
mission3 = [137.6,  72.9, 67.9, 46.1];
mission4 = [117.8, 103.4, 79.1, 89.4];

comms(1,:,:) = [0,0; 1, 6; 0,  5; 0,0];
comms(2,:,:) = [0,0; 0, 2; 5, 14; 0,0];
comms(3,:,:) = [0,0; 0, 3; 4,  3; 0,0];
comms(4,:,:) = [0,0; 1, 4; 0,  6; 0,0];


vals = [ mission1; mission2; mission3; mission4];

names = ['Mission 1'; 'Mission 2'; 'Mission 3'; 'Mission 4'];

figure
subplot(2,1,1)
Graph( 'Field Trial', names, vals, {'Dead Reckoning'; 'TBN'; 'dec-TBN Full'; 'dec-TBN Planned'});
subplot(2,1,2)
plotBarStackGroups(comms, names)




%% Simulation 4 AUVs
clear all
close all
clc

vals = [8.54, 70.3, 11.66, 17.44, 14.82, 12.49, 10.59, 9.31, 8.93, 9.39, 10.15];

names = {'Policy', 'Dead Reckoning', 'TBN', 'Full (100%)', '80%', '60%', '40%', '20%', '10%', '5%', '2.5%'};

Graph( '4 AUV Simulation', names, vals)




%% Simulation 2 AUVs
clear all
close all
clc

vals = [10.53, 70.4, 12.36, 18.85, 16.89, 15.81, 14.44, 12.30, 11.52, 11.17, 11.27];

names = {'Policy', 'Dead Reckoning', 'TBN', 'Full (100%)', '80%', '60%', '40%', '20%', '10%', '5%', '2.5%'};

Graph( '2 AUV Simulation', names, vals)





%% Do plotting
function Graph( title_, xNames, vals, lNames)


%figure('name', sprintf('%s Results', title_), 'numbertitle', 'off')

b = bar(vals);

for ii = 1: size(vals,2)
    xtips1 = b(ii).XEndPoints;
    ytips1 = b(ii).YEndPoints;
    labels1 = string(b(ii).YData);
    text(xtips1,ytips1,labels1, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', ...
        'fontsize', 14)
end

ax = gca;
ax.FontSize = 16;               % Set font size first or else you will loose the rest of the formatting
ax.XTickLabel = xNames;          % Use mission names
% ax.XTickLabelRotation = 45;

if nargin == 4
    l = legend(lNames);
    l.FontSize = 16;
    l.Location = 'northeastoutside';
end

title( sprintf("%s: Average Error", title_), 'fontsize',26)
ylabel('Average Error [m]', 'fontsize', 20)


end