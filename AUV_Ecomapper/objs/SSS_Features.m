% Side Scan Sonar Features
% Jake Anderson
% 10/1/2019

classdef SSS_Features
    
    properties
        Features
        MatchingData
        Header
        
    end
    
    methods (Static)
        
        function obj = SSS_Features, warning on backtrace
            
            obj.Features.ImagePatch(imageSize, imageSize) = 0;
            obj.Features.ImageName = 'Name Goes Here';
            obj.Features.Centroid  = [0,0];
            obj.Features.UTMzone   = "";
            obj.Features.Boundary  = [];
            obj.Features.Cov       = zeros(2,2);
            obj.Features.Radius    = 0;
            obj.Features.Index     = [0,0];
            obj.Features.Timestamp = datetime;
            
            obj.MatchingData.Probs = [];
            obj.MatchingData.Truth = false;
            
            obj.Header = [ "ImagePatch: An image patch centered on the feature";
                           " ImageName: Name of the saved image patch jpeg";
                           "  Centroid: UTMs of the feature's locaiton";
                           "  Boundary: Ploygon that outlines the feature in the master Image";
                           "       Cov: Covariance of the Feature's locaiton relative to the vehicle";
                           "    Radius: Feature's distance away from the vehicle; + is satboard, - is port";
                           "     Index: Center of the feature on the master image";
                           " TimeStamp: Time stamp of the feature's Index";
                           "     Prods: Feature matching probabilites";
                           "     Truth: Indicates ture feature matches"];
            
        end
        
        
    end
    
    
    methods
        
        % Feature Position in lat - lon
        function [lat, lon] = Feature_latlon(obj)
            
            features = cat(1,obj.Features.Centroid);
            
            if ~isempty(features)
                utmZone = obj.Features.UTMzone;
                [lat, lon] = utm2deg(features(:,1), features(:,2), repmat(utmZone, size(features,1),1));
            else
                lat = [];
                lon = [];
            end
            
        end
        

    end
    
end