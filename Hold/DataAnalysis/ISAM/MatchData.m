%% Get Feature Matching Data
%  Jacob Anderson
%  8/7/2019

classdef MatchData
   
    properties
        TureMatches
        Names
        Probs
        Truth
        XY
        Err
        Indx
        Time
    end
    
    properties (Access = private)

    end
    
    methods
        
        % Instanciate Class by Get Matching Predictions -------------------------------------------------
        function obj = MatchData(ui)
            
            
            file = ui.getPath('FeatureMatchingData','txt','file');          % Get the name and file path of preditions txt file
            
            fileID = fopen(file,'r');
            matchData = textscan(fileID, '%s %s %f', 'Delimiter',',');
            fclose(fileID);
            
            names = unique([matchData{1}; matchData{2}]);
            
            Prob = zeros(size(names,1));
            
            names1 = matchData{1};
            names2 = matchData{2};
            prob   = matchData{3};
            
            for ii = 1 : size(matchData{1})
                
                row = find(strcmp(names1(ii), names));
                col = find(strcmp(names2(ii), names));
                
                Prob(row, col) = prob(ii);
                Prob(col, row) = prob(ii);
                
            end
            
            obj.Probs = Prob;
            obj.Names = names;
            obj.Truth = false(size(Prob));
            obj.XY    = NaN(size(names,1),2);
            obj.Err   = NaN(size(names,1),2);
            obj.Indx  = NaN(size(names,1),1);
            obj.Time  = NaT(size(names,1),1);
            
        end
        
        
        % Get True Matches ------------------------------------------------------------------------------
        function obj = GetTrueMatches(obj, starFeatures, portFeatures, dist_threshold)
            
            fprintf("\n\nFinding True Feature Matches\nAllowing a maximum seperation of: %f\n\n", dist_threshold)
            
            obj.Truth = false(size(obj.Probs));                             % Re-set the Truth array
            
            centroid   = [cat(1, portFeatures.Centroid);  cat(1, starFeatures.Centroid) ];
            imagePatch = [{portFeatures.imagePatch}';      {starFeatures.imagePatch}'     ];
            imageName  = [cat(1, portFeatures.imageName); cat(1, starFeatures.imageName)];
            timestamp  = [cat(1, portFeatures.timestamp); cat(1, starFeatures.timestamp)];
            
            indexes = 1: size(centroid,1);
            
            
            jj = 1;
            for ii = 1: size(centroid,1)
                
                ind = false(size(centroid,1),1);
                
                ind(ii) = true;
                
                currentFeat = centroid(ind,:);
                otherFeat   = centroid(~ind,:);
                
                index = indexes(~ind);
                
                dist = vecnorm( (otherFeat - currentFeat)' );
                
                match = dist <= dist_threshold;
                
                if any(match)
                    
                    inds = index(match);
                    
                    for aa = 1: numel(inds)
                        
                        [row, col] = obj.GetIndex( imageName( ind,:), imageName( inds(aa),:) );
                        
                        obj.Truth(row, col) = true;
                        obj.Truth(col, row) = true;
                        
                        obj.TureMatches.A(jj).Centroid   = centroid(  ind,:);
                        obj.TureMatches.A(jj).imageName  = imageName( ind,:);
                        obj.TureMatches.A(jj).imagePatch = imagePatch{ind,:};
                        obj.TureMatches.A(jj).timestamp  = timestamp( ind,:);
                        
                        obj.TureMatches.B(jj).Centroid   = centroid(  inds(aa),:);
                        obj.TureMatches.B(jj).imageName  = imageName( inds(aa),:);
                        obj.TureMatches.B(jj).imagePatch = imagePatch{inds(aa),:};
                        obj.TureMatches.B(jj).timestamp  = timestamp( inds(aa),:);
                        
                        jj = jj+1;
                        fprintf('%s <--> %s\n', imageName(ind,:), imageName(inds(aa),:))
                    end
                end
                
                
            end
            
        end
        
        
        % Assign location to Feature --------------------------------------------------------------------
        function obj = AddLocation(obj, features, err)
            
            names    = cat(1,features.names); 
            location = cat(1,features.xy);
            times    = cat(1, features.times);
            
            for ii = 1: size(names,1)
                ind = strcmp(names(ii,:), obj.Names);
                obj.XY(ind,:)  = location(ii,:);
                obj.Err(ind,:) = err(ii);
                obj.Time(ind) = times(ii);
            end
            
        end
        
        
        % Assign Index ----------------------------------------------------------------------------------
        function obj = AddIndex(obj, names, indx)
            
            for ii = 1: size(names,1)
                ind = strcmp(names(ii,:), obj.Names);
                obj.Indx(ind) = indx + ii;
            end
            
        end
            
        
        % Get Potential Matches -------------------------------------------------------------------------
        function matches = PotnetialMatch(obj, features, prob_threshold)
            
            names   = cat(1,features.names);
            feat_xy = cat(1, features.xy);
            cov_    = features.cov;
            
            self = ismember(obj.Names, names);                              % The the current feature's index
            ind  = ~isnan(obj.Indx(:));                                     % Get the indexes of features that have a valid location
            ind(self) = false;                                              % Do not compare the feature to itself
            
            if sum(ind) < 2                                                 % There isn't anything to compare to
                matches = [];
                return
            else
                matches(size(names,1)).Feature = '';                        % Pre-alocate memeory
                matches(size(names,1)).Matches = '';
                matches(size(names,1)).Indx    = [];
                matches(size(names,1)).Valid   = false;
            end
            
            
            count = 1;
            for ii = 1: size(names,1)                                       % Iterate throught the new features
                
                self = strcmp(names(ii,:), obj.Names);                      % The the current feature's index
                
                self_xy    = obj.XY(self,:);                                % Current Feature's X Y location
                self_err   = obj.Err(self,:);                               % Error in location
                self_time  = obj.Time(self,:); 
                
                comp_xy    = obj.XY(ind,:);                                 % Data corsponding to the potentila matches
                comp_err   = obj.Err(ind,:);
                comp_names = obj.Names(ind,:);
                comp_inds  = obj.Indx(ind,:);
                comp_time  = obj.Time(ind,:);
                
                dist = vecnorm((comp_xy - self_xy)');                       % Distances to the other features
                err  = vecnorm((comp_err + self_err)');
                
                in = dist <= err & self_time' > comp_time';                 % Distance thresholde
                
                comp_xy    = comp_xy(in,:);
                comp_names = comp_names(in,:);                              % Reduce data
                comp_inds  = comp_inds(in,:);
                
                prob  = obj.GetProb( names(ii,:), comp_names );             % Get probabilities of the proposed matches
                truth = obj.GetTruth(names(ii,:), comp_names );
                
                in = prob >= prob_threshold;                                % Probabiltiy threshold
 
                if any(in)                                                  % We think there is a match
                    matches(count).Feature = names(ii,:);                   % Assign data
                    matches(count).Matches = comp_names{in,:};
                    matches(count).XY      = comp_xy(in,:);
                    matches(count).XY_feat = feat_xy(ii,:);
                    matches(count).Indx    = comp_inds(in);
                    matches(count).Valid   = truth(in);
                    matches(count).cov     = cov_(:,:,ii);
                    count = count+1;
                end
            end
            
            matches(count:end) = [];                                        % Eliminate extra pre-alocated Space
                       
        end
        
        
        % Reset Feature matching Data for new simulation ------------------------------------------------
        function obj = Reset(obj)
            obj.XY    = NaN(size(obj.Names,1),2);
            obj.Err   = NaN(size(obj.Names,1),2);
            obj.Indx  = NaN(size(obj.Names,1),1);
            obj.Time  = NaT(size(obj.Names,1),1);
        end
        
        
        % Display Image Patches from the true Matches
        function DisplayMatches(obj)
            
            for ii = 1:size(obj.TureMatches.A,2)
                
                fig_name = sprintf("Match %d",ii);
                figure('name', fig_name, 'numbertitle', 'off')
                subplot(1,2,1)
                imshow(obj.TureMatches.A(ii).imagePatch)
                title(obj.TureMatches.A(ii).imageName)
                
                subplot(1,2,2)
                imshow(obj.TureMatches.B(ii).imagePatch)
                title(obj.TureMatches.B(ii).imageName)
            end

        end
        
        
        % Posible Matches by probabiltiy threshold
        function MatchByProb(obj, threshold, AllFeatures)
            
            [row,col] = find(obj.Probs >= threshold);
            
            fprintf("\n\nPosible Matches based on a %f threshold\n",threshold)
            fprintf("%s <--> %s\n",obj.Names{row}, obj.Names{col})
            
            if nargin < 3, return, end
                       
            in_1 = ismember(AllFeatures.imageName, cat(1, obj.Names{row}), 'rows');
            in_2 = ismember(AllFeatures.imageName, cat(1, obj.Names{col}), 'rows');
            
            matchImage1 = AllFeatures.imagePatch(in_1);
            matchImage2 = AllFeatures.imagePatch(in_2);
            
            s = size(matchImage1,1);
            
            for ii = 1:s
                
                figure('Name',"Proposed Image matches")
                
                subplot(1,2,1)
                imshow(matchImage1{ii})
                
                subplot(1,2, 2)
                imshow(matchImage2{ii})
                
            end
            
        end
        
    end
    
    
    
    methods (Access = private)
        
                % Get matching Probability ----------------------------------------------------------------------
        function prob = GetProb(obj, name1, name2 )
            
            if isempty(name1) || isempty(name2)
                prob = 0;
                return
            end
            
            prob = zeros(size(name1,1),size(name2,1));
            
            for ii = 1: size(name1,1)
                for jj= 1: size(name2,1)
                    
                    [row, col] = obj.GetIndex(name1(ii,:), name2(jj,:));
                    prob(ii,jj) = obj.Probs(row,col);
 
                end
            end
            
        end
        
        
        % Get Ground Truth Matching data ----------------------------------------------------------------
        function truth = GetTruth(obj, name1, name2 )
            
            if isempty(name1) || isempty(name2)
                truth = false;
                return
            end
            
            truth = zeros(size(name1,1),size(name2,1));
            
            for ii = 1: size(name1,1)
                for jj= 1: size(name2,1)
                    
                    [row, col] = obj.GetIndex(name1(ii,:), name2(jj,:));
                    truth(ii,jj) = obj.Truth(row,col);
                  
                end
            end
           
        end
    
        
        % Get index of probabilty matrix ------------------------------------------------------------------
        function [row, col] = GetIndex(obj, name1, name2)
            row = find(strcmp(name1, obj.Names));
            col = find(strcmp(name2, obj.Names));
        end
        
    end
end
