%% Manifest of features being sent
classdef Manifest
    
    properties (Access = private)
        names
        values
        error
        centroids
    end
    
    methods      
        
        % Constructor
        function obj = Manifest( starFeatures, portFeatures)
            
            starNames = {starFeatures.imageName};
            portNames = {portFeatures.imageName};
            
            keys = [starNames'; portNames'];
            
            obj.names     = keys;
            obj.values    = cell(size(keys));
            obj.error     = zeros( size(keys,1), 2);
            obj.centroids = zeros( size(keys,1), 2);
            
        end
        
        
        % Update the entries in the Manifest
        function obj = UpdateManifest(obj, features)
            
            names_   = features.names; 
            location = features.xy; 
            count    = features.count;
            err      = features.std;
            
            for ii = 1: size(names_,1)
                
                indx = ismember(obj.names, names_(ii,:));
                
                if isempty(obj.values{indx})
                    obj.values{indx} = count + ii;
                    obj.error(indx,:) = err;
                    obj.centroids(indx,:) = location(ii,:);
                    
                else
                    fprintf('\n\nManifest indicates multiple entries for item %s\n\n', names_)
                    obj.values{indx} = [obj.values{indx}, count + ii];
                end
            end
            
        end
        
        
        % Retrive entries from the manifest
        function [ids, error] = GetID(obj, names )
           
            ids = cell(size(names));
            err = zeros(size(names,1),2);
            
            for ii = 1: size(names,1)
                indx = ismember(obj.names, names(ii));
                ids(ii) = obj.values(indx);
                err(ii,:) = obj.error(indx,:);
            end
            
            if nargout > 1
                error = err;
            end
            
        end
        
        
        % Get near by Entries
        function nearby = Neighbors(obj, features, sig)
            
            nearby.Names = '';
            nearby.Ids   = {};
            nearby.Err   = [];
            
            jj = 1;
            for ii = 1: size(features.xy,1)
                dist = vecnorm( (obj.centroids - features.xy(ii,:))') + vecnorm(obj.error');
                indx = dist <= sig * norm(features.std);
                
                disp("Threhold")
                disp(sig * norm(features.std));
                
                disp("Min Dist")
                disp(min(dist))
                
                if any(indx)
                    nearby.Names(jj,:) = obj.names{indx};
                    nearby.Ids(jj,:)   = obj.values(indx);
                    nearby.Err(jj,:)   = obj.error(indx,:);
                    jj = jj+1;
                end
            end
            
        end
        
    end
end

