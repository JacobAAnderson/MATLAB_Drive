% Model Acoustic Modems
% Jacob Anderson
% RDML, OSU, Corvallis OR.
% Sept 10, 2020

close all
clear all
clc


%% Get Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                      % Indicate whether new data should be selescted by the user


filters = {'GPS_fix', 0};                                                   % Data filter {"Data_Field", value to be filterred out}
%           'ALT_ALTITUDE', 0};        

[RT, geotiff, bathy] = GetRiptides(ui, 2, {'210','216'}, filters);          % Choose data files with gui and return data structure

RT_ = RT;                                                                   % Keep a clean coppy of the vehicle data


%% Prep Data
clc
% missions = [3,2];         % Test Lake Mission Set
missions = [3,4];           % Fosters Lake Mission Set

for v = 1: numel(RT)
    
    RT(v) = RT(v).Select_Mission( missions(v) );                            % Choose Mission
    RT(v) = RT(v).Model_VehicleSpeed;                                       % Calculate corrected speed
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths
    RT(v) = RT(v).Add_ParticleFilter(bathy);                                % Prep Particle filters
    
    RT(v).Disp_Manifest;                                                    % Display Vehicle Manifest
    
end

clear missions v




%% ---- Model the latency in acoustic communications ---
[mu, sig] = Model_AcommsLatency( RT );



%% --- Show Acoustic Communications ---
[~] = Plot_AcousticCommunications(RT, geotiff);



