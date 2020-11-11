%% Run Stats on PF_SLAM results
close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

file = ui.GetFile('RMS_Data','mat');                                        % GUI to get the name and file path of a file

if isempty(file), return, end                                               % Check if the Data Input box was cnaceled

load(file)

clear file

%results.ClearEmpty;
%%
% results.Plot_RMS;

s = size(results.RMS_err(1).PF, 1);

gtrms = rmse(results.GT(1:s,:) - results.DR(1:s,:));

pf = cat(2, results.RMS_err.PF);
pfave = mean(pf,2);


pfslam = cat(2, results.RMS_err.PF_SLAM);
pfslamave = mean(pfslam,2);

figure
plot(gtrms,'r')
hold on
plot(pfave,'b')
plot(pfslamave,'c')
hold off
ylabel("RMSE [meters]")
xlabel("Time Steps")
legend("Dead Reckoning", "Particle Filter", "SLAM", 'Location','northwest')
title("Average RMSE over 50 trials")



%% Matching
clear dr bathyPercent count iter
close all
clc

file = ui.GetFile('ISAM_Data','mat');                                       % GUI to get the name and file path of a file
if isempty(file), return, end                                               % Check if the Data Input box was cnaceled

load(file)

clear file newData

%%
clc
close all
clear tureMatches

AllFeatures.Centroid   = [cat(1, portFeatures.Centroid);  cat(1, starFeatures.Centroid) ]; 
AllFeatures.imagePatch = [{portFeatures.imagePatch}';      {starFeatures.imagePatch}'     ]; 
AllFeatures.imageName  = [cat(1, portFeatures.imageName); cat(1, starFeatures.imageName)]; 
AllFeatures.timestamp  = [cat(1, portFeatures.timestamp); cat(1, starFeatures.timestamp)]; 


matchData.MatchByProb(0.6, AllFeatures);

indexes = 1: size(AllFeatures.Centroid,1);


figure('name', 'All Features', 'numbertitle', 'off')
plot(AllFeatures.Centroid(:,1), AllFeatures.Centroid(:,2), '+k')
hold on
plot( match1(:,1), match1(:,2), '+b')
plot( match2(:,1), match2(:,2), '+r')


jj = 1;
for ii = 1: size(AllFeatures.Centroid,1)
    
    ind = false(size(AllFeatures.Centroid,1),1);
    
    ind(ii) = true;
    
    currentFeat = AllFeatures.Centroid(ind,:);
    otherFeat   = AllFeatures.Centroid(~ind,:);
    
    index = indexes(~ind);
    
    dist = vecnorm( (otherFeat - currentFeat)' );
    
    match = dist <= 2;
    
    if any(match)
        
        inds = index(match);
        
        tureMatches.A(jj).Centroid   = AllFeatures.Centroid(  ind,:);
        tureMatches.A(jj).imageName  = AllFeatures.imageName( ind,:);
        tureMatches.A(jj).imagePatch = AllFeatures.imagePatch(ind,:);
        tureMatches.A(jj).timestamp  = AllFeatures.timestamp( ind,:);
        
        tureMatches.B(jj).Centroid   = AllFeatures.Centroid(  inds,:);
        tureMatches.B(jj).imageName  = AllFeatures.imageName( inds,:);
        tureMatches.B(jj).imagePatch = AllFeatures.imagePatch(inds,:);
        tureMatches.B(jj).timestamp  = AllFeatures.timestamp( inds,:);
        
        
        jj = jj+1;
        fprintf('%s <--> %s\n', AllFeatures.imageName(ind,:), AllFeatures.imageName(inds,:)')
        
%         plot(AllFeatures.Centroid(ind,1), AllFeatures.Centroid(ind,2) ,'o')
        
    end
    

end

match1 = matches(1).featureLocation;
match2 = matches(2).featureLocation;


for ii = 1: size(match1,1) 
    plot( [match1(ii,1), match2(ii,1)], [match1(ii,2), match2(ii,2)], '+')
end

hold off

clear ii jj match dist ind inds indx index indexes currentFeat otherFeat


%% 
clc

matches(1).valid = false(size(match1,1),1);

% Find Positive matches
for ii = 1: size(match1,1)
    
    dist = vecnorm( match1(ii,:) - match2(ii,:) );
    
    if dist < 2
        matches.valid(ii) = true;
    end

end

