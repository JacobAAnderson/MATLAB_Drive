% Multi-Robot AUV
% Jacob Anderson
% 1/30/2020

function auv = MR_AUV1( name, bathy, wpFile_self, varargin)

auv = AUV1( name, bathy, wpFile_self);

% ___ Communications Planning _____________________________________________
auv.Planning = CommunicationPlanning(name);

pf = auv.Nav.PF;
pf = pf.SetNoise('speed',  'Normal', 0, 1.0);                               % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('compass','Normal', 0, 10);                                % Set Comapss Noise
pf = pf.SetNoise('alt',    'Normal', 0, 0.5);                               % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('acoms',  'Normal', 0, 0.5);                               % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )


rad = @(x) x*pi/180;

speedNoise   = makedist('Normal','mu',0.2,'sigma',0.5);
compassNoise = makedist('Normal','mu',rad(3),'sigma',rad(10));

% Simulate Self
path = csvread(wpFile_self);
auv.Planning = auv.Planning.Add_Vehicle(pf, path, 'self', 3, 3, speedNoise, compassNoise);

% Simulate Other Vehicles
for ii = 1: numel(varargin)
    path = csvread(varargin{ii});
    pf = pf.SetLocation(path(1,1:2), eye(2) * 3);
    auv.Planning = auv.Planning.Add_Vehicle(pf, path, sprintf('vehicle%d',ii), 3, 3, speedNoise, compassNoise);
end

end

