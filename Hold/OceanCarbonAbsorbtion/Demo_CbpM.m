% Demo CbPm model
% Jacob Anderson
% 2/18/2020

close all
clear all
clc

fprintf('\nLoading Data\n')
%% Load Data

od = OceanData;

if exist('/home/jake/Desktop/A20031522017181_par_9km.nc', 'file')
    
    od = od.AddData('par',          '/home/jake/Desktop/A20031522017181_par_9km.nc');
    od = od.AddData('Kd_490',       '/home/jake/Desktop/A20031522017181_Kd490_9km.nc');
    od = od.AddData('chl_gsm',      '/home/jake/Desktop/A20031522017181_GSM_chl_gsm_9km.nc');
    od = od.AddData('bbp_443_gsm',  '/home/jake/Desktop/A20031522017181.L3m_MC_GSM_bbp_443_gsm_9km.nc');
    od = od.AddData('mld',          '/home/jake/Desktop/mld_2003152.hdf');
    od = od.AddData('zno3',         '/home/jake/MATLAB-Drive/OceanCarbonAbsorbtion/src/monthly_nitracline.mat');
    

    
elseif exist('/Users/jake/Desktop/A20031522017181.L3m_MC_PAR_par_9km.nc','file')
    
    od = od.AddData('par',          '/Users/jake/Desktop/A20031522017181.L3m_MC_PAR_par_9km.nc');
    od = od.AddData('Kd_490',       '/Users/jake/Desktop/A20031522017181.L3m_MC_KD490_Kd_490_9km.nc');
    od = od.AddData('chl_gsm',      '/Users/jake/Desktop/A20031522017181.L3m_MC_GSM_chl_gsm_9km.nc');
    od = od.AddData('bbp_443_gsm',  '/Users/jake/Desktop/A20031522017181.L3m_MC_GSM_bbp_443_gsm_9km.nc');
    od = od.AddData('mld',          '/Users/jake/Desktop/mld.2003152.hdf');
    od = od.AddData('zno3',         '/Users/jake/MATLAB-Drive/OceanCarbonAbsorbtion/src/monthly_nitracline.mat');
    
    
    
elseif exist('/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/par/A20031522017181.L3m_MC_PAR_par_9km.nc', 'file')
    % Server Stuff
    % addpath(genpath('/data3/Jacob/OceanCarbonAbsorbtion'))
    %
%     A = exist('/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/par/A20031522017181.L3m_MC_PAR_par_9km.nc', 'file');
%     A = exist('/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/Kd/A20031522017181.L3m_MC_KD490_Kd_490_9km.nc', 'file');
%     A = exist('/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/chl_gsm/A20031522017181.L3m_MC_GSM_chl_gsm_9km.nc', 'file');
%     A = exist('/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/bbp/A20031522017181.L3m_MC_GSM_bbp_443_gsm_9km.nc', 'file');
%     A = exist('/data1/mld/mld030.hycom/monthly/mld.2003152.hdf', 'file');
%     A = exist('/data3/Jacob/OceanCarbonAbsorbtion/src/monthly_nitracline.mat', 'file');
    
    
    od = od.AddData('par',          '/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/par/A20031522017181.L3m_MC_PAR_par_9km.nc');
    od = od.AddData('Kd_490',       '/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/Kd/A20031522017181.L3m_MC_KD490_Kd_490_9km.nc');
    od = od.AddData('chl_gsm',      '/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/chl_gsm/A20031522017181.L3m_MC_GSM_chl_gsm_9km.nc');
    od = od.AddData('bbp_443_gsm',  '/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/bbp/A20031522017181.L3m_MC_GSM_bbp_443_gsm_9km.nc');
    od = od.AddData('mld',          '/data1/mld/mld030.hycom/monthly/mld.2003152.hdf');
    od = od.AddData('zno3',         '/data3/Jacob/OceanCarbonAbsorbtion/src/monthly_nitracline.mat');
    
else
    disp('Where am I??????')
end

%% Modify Data
fprintf('\nModifying The Data\n')

od = od.Modify('mld',  @(x) flipud(x));                                     % Flip mld data to match the other data sets
od = od.Modify('zno3', @(x) x(:,:,6));                                      % Select month 
od = od.Modify('zno3', @(x) flipud(x));                                     % Flip zno3 data to match the other data sets

od = od.World_Coordinates({'par','Kd_490','chl_gsm','bbp_443_gsm'});                     % Create Unified World Coordinate System
od = od.Fit2World({'mld', 'zno3'});                                         % Fit the rest of the data parmaeters to the lat lon grid

od = od.Daylength('mld');                                                   % Calculate day length for world coordinate frame


%% Plot some of the data
% fig1 = od.Plot('par');
% fig2 = od.Plot('mld');
% fig3 = od.Plot('zno3');


%% Calculate Daylighe
fprintf('\nProcessing Model\n')

% flds = fields(od.Data.bbp_443_qaa);
% disp(flds)

scale = 0.25;
% Assigne Values
par  = od.Value('par',          scale);                                                     
k490 = od.Value('Kd_490',       scale);
mld  = od.Value('mld',          scale);
zno3 = od.Value('zno3',         scale);
dl   = od.Value('dl',           scale);
chl  = od.Value('chl_gsm',      scale);
bbp  = od.Value('bbp_443_gsm',  scale);


ppm = PPM(bbp, chl, par, k490, mdl, dl);

% [ppz, mu, chlz, Cz, PARz, Kz, Zeu,Ig]=CbPM_fxn(par,K490,mld,zno3,dl,chl,bbp);


%% One at a time -- Maybee???

% %matObj = matfile("Results/CbPM_fxn.mat"); 
% 
% 
% m = matfile('Results_CbPM_fxn.mat','Writable',true);
% 
% % Prealocate Memory
% 
% [r,c] = size(par);
% 
% m.ppz(r,c,300)  = 0;
% m.mu(r,c,300)   = 0;
% m.chlz(r,c,300) = 0;
% m.Cz(r,c,300)   = 0;
% m.PARz(r,c,300) = 0;
% m.Kz(r,c,300)   = 0;
% m.Zeu(r,c,300)  = 0;
% m.Ig(r,c,300)   = 0;
% 
% %save('Results_CbPM_fxn.mat', 'ppz', 'mu', 'chlz', 'Cz', 'PARz', 'Kz', 'Ig')
% 
% for r = 1: size(par,1)
%     for c = 1:size(par,2)
% 
%         disp([r,c])
%         
%         try
%         
%         [ppz_, mu_, chlz_, Cz_, PARz_, Kz_, Zeu_, Ig_] = CbPM_fxn( par(r,c), K490(r,c), mld(r,c), zno3(r,c), dl(c), chl(r,c), bbp(r,c));
%     
%         m.ppz(r,c,:)  = shiftdim( real(ppz_),  -1);
%         m.mu(r,c,:)   = shiftdim( real(mu_),   -1);
%         m.chlz(r,c,:) = shiftdim( real(chlz_), -1);
%         m.Cz(r,c,:)   = shiftdim( real(Cz_),   -1);
%         m.PARz(r,c,:) = shiftdim( real(PARz_), -1);
% %        m.Kz(r,c,:)   = shiftdim( real(Kz_),   -1);
% %         m.Zeu(r,c,:)  = Zeu_;
%         m.Ig(r,c,:)   = shiftdim( real(Ig_), -1);
% 
%         catch
%             
%         end
%     
%     end
% end
% 
% % save('Results_CbPM_fxn.mat', 'ppz', 'mu', 'chlz', 'Cz', 'PARz', 'Kz', 'Ig') 





