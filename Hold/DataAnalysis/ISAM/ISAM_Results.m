% ISAM Results
% Jake Anderson
% 10/3/2019

classdef ISAM_Results
    
    properties
        GT
        DR
        Paths
        RMS_err
        MapDim
        Param
        
    end
    
    methods
        
        function obj = ISAM_Results(num)
            
            obj.Paths(num).PF      = [];
            obj.Paths(num).SLAM    = [];
            obj.Paths(num).PF_SLAM = [];
                        
            obj.RMS_err(num).PF      = [];
            obj.RMS_err(num).SLAM    = [];
            obj.RMS_err(num).PF_SLAM = [];
            obj.Param(num) = 0;
            
        end
        
        function obj = AddGroundTruth(obj, gt)
            obj.GT = gt;
        end
        
        function obj = AddDearReconing(obj, dr) 
            obj.DR = dr;
        end
        
        function obj = AddMapDim(obj, size_, range)
            obj.MapDim.size = size_;
            obj.MapDim.range = range;
        end
        
        function obj = Truncate(obj,ind)
            obj.Paths(ind:end)   = [];
            obj.RMS_err(ind:end) = [];
            obj.Param(ind:end)   = [];
        end
        
        
        function obj = ClearEmpty(obj)
            
            out = arrayfun(@isempty, [obj.RMS_err.PF]);
            ind = size(out,2) + 1;
            
            obj.Paths(ind:end)   = [];
            obj.RMS_err(ind:end) = [];
            obj.Param(ind:end)   = [];
        end
        
        
        function Plot_RMS(obj)
            
            % --- Average out eh RMS Error at the various resolutions over the different trials -------
            
%             bathyPercent = unique( cat(1, obj.Param) );                   % Get all of the bathymetry resolutions
            bathyPercent = 0.1: 0.05: 0.4;
            bathyPercent(bathyPercent == 0) = [];                           % Get Rid of 0
            
            bathyRes = vecnorm( (obj.MapDim.range' ./ ([bathyPercent;bathyPercent] .* fliplr( obj.MapDim.size)')) );
            
            s = max(size(bathyRes));
            
            rmsErr(s).pf     = zeros(size(obj.RMS_err(1).PF));              % Pre-alocate memory
            rmsErr(s).pf_ave = 0;
            rmsErr(s).pf_max = 0;
            
            rmsErr(s).slam     = zeros(size(obj.RMS_err(1).PF));
            rmsErr(s).slam_ave = 0;
            rmsErr(s).slam_max = 0;
            
            rmsErr(s).pf_slam     = zeros(size(obj.RMS_err(1).PF));
            rmsErr(s).pf_slam_ave = 0;
            rmsErr(s).pf_slam_max = 0;
            
            rmsErr(s).bathy_percent = 0;
            rmsErr(s).bathy_res = 0;
            
            rmsErr(1).batyh_range = obj.MapDim.range;
            rmsErr(1).batyh_size  = obj.MapDim.size;
            
            
            for ii = 1: s
                
                ind = ii:s:size(obj.RMS_err,2);
                
                rmsErr(ii).bathy_res = bathyRes(ii);
                rmsErr(ii).bathy_percent = bathyPercent(ii);
                
                rmsErr(ii).pf      = mean(cat(2,obj.RMS_err(ind).PF),2);
                rmsErr(ii).slam    = mean(cat(2,obj.RMS_err(ind).SLAM),2);
                rmsErr(ii).pf_slam = mean(cat(2,obj.RMS_err(ind).PF_SLAM),2);
                
                rmsErr(ii).pf_ave = mean(rmsErr(ii).pf);
                rmsErr(ii).pf_max =  max(rmsErr(ii).pf);
                
                rmsErr(ii).slam_ave = mean(rmsErr(ii).slam);
                rmsErr(ii).slam_max =  max(rmsErr(ii).slam);
                
                rmsErr(ii).pf_slam_ave = mean(rmsErr(ii).pf_slam);
                rmsErr(ii).pf_slam_max =  max(rmsErr(ii).pf_slam);
                
            end
            
            % --- Get deadreckoning RMS -----------------------------------
            dr.rms = rms(obj.DR -obj.GT);
            dr.rms_ave = mean( dr.rms );
            dr.rms_ave = repmat(dr.rms_ave, size(bathyPercent));
            
            dr.rms_max = max( dr.rms );
            dr.rms_max = repmat(dr.rms_max, size(bathyPercent));
            
            
            % --- Plottting -----------------------------------------------
            figure('name','Ave RMS Error Stats', 'numbertitle', 'off');
            plot( bathyPercent, cat(1, rmsErr.pf_ave),'b');
            hold on
            plot( bathyPercent, cat(1, rmsErr.slam_ave),'m');
            plot( bathyPercent, cat(1, rmsErr.pf_slam_ave),'g');
            plot( bathyPercent, dr.rms_ave,'r');
            hold off
            title('Average RMS Path Error')
            xlabel('Bathymetery Map Resolution [% of master Map]')
            ylabel('RMS Error [m]')
            legend('PF', 'SLAM', 'PF-SlAM', 'DR')
            
            
            figure('name','Max RMS Error Stats', 'numbertitle', 'off')
            plot( bathyPercent, cat(1, rmsErr.pf_max),'b');
            hold on
            plot( bathyPercent, cat(1, rmsErr.slam_max),'m');
            plot( bathyPercent, cat(1, rmsErr.pf_slam_max),'g');
            plot( bathyPercent, dr.rms_max,'r');
            hold off
            title('Max RMS Path Error')
            xlabel('Bathymetery Map Grid size [% of master Map]')
            ylabel('RMS Error [m]')
            legend('PF', 'SLAM', 'PF-SlAM', 'DR')
            
        end
        
    end
end