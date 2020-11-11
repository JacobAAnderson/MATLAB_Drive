

function asv = ASV1(bathymetry, start)

asv = R0b0t('ASV');

asv = asv.SetUncertainty('Speed',   'Normal', 0, 0.25);
asv = asv.SetUncertainty('Compass','Normal', 0, 10);

if nargin < 2
    asv = asv.SetLocation(362000, 3701500);
else
    asv = asv.SetLocation(start(1), start(2));
end

bathymetry = bathymetry.ReduceResolution(0.25);

asv = asv.Add_Map(bathymetry.AsMapp);

asv = asv.Add_Sensor("Bathy", 'Normal', 0, 0.01);
asv = asv.Add_Sensor("Salt",  'Normal', 0, 0.01);
asv = asv.Add_Sensor("GPS",   'Normal', 0, 1.5);

end