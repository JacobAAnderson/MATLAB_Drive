%% Water Feature Model

function wf = WaterFeature3(bathymetry, plotIO)

XX = bathymetry.Easting;
YY = bathymetry.Northing;
ZZ = bathymetry.Elevation;

wf =  WaterFeature(XX, YY, ZZ, "Sim Feature");


% First Feature       Mu         Sigma   Rho   S
wf = wf.Add_Gaussian([-2,   3], [4, 3],  0.5,  5);                          % First  Gaussian
wf = wf.Add_Gaussian([ 2,   4], [4, 2], -0.4,  3);                          % Second Gaussian
wf = wf.Add_Gaussian([-1.5, 7], [2, 4], -0.5,  3);                          % Third  Gaussian
wf = wf.Add_Gaussian([-2,   3], [2, 2], -0.0,  2);                          % Fourth Gaussian

% % Second Feature
% wf = wf.Add_Gaussian([-3, -3], [1.3, 1.0],  0.0, 1.3);                      % First  Gaussian
% wf = wf.Add_Gaussian([-1, -4], [1.0, 1.3],  0.0, 1.3);                      % Second Gaussian


% Tird Feature
wf = wf.Add_Gaussian([ 5, -2], [1.3, 1.0],  0.0,  1);                      % First  Gaussian
wf = wf.Add_Gaussian([ 5, -3], [4.0, 5.3], -0.5, 20.0);                    % Second Gaussian



wf = wf.MakeMap;


% wf = wf.SetUncertainty('Speed',    'Normal',0.0, 0.0);
% wf = wf.SetUncertainty('Compass',  'Normal',0.0, 0.0);
% wf = wf.SetUncertainty('Spreading','Normal',0.0, 0.0);

if nargin > 1 && plotIO, wf.PlotMap; end

end

