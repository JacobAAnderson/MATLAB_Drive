

function map = LoadBathymetry_asMap(bathFile)

bathymetry = LoadBathymetry_asMap(bathFile);                                      % Load Bathymetry Map

XX = bathymetry.Easting;
YY = bathymetry.Northing;
ZZ = bathymetry.Elevation;

map = Mapp(XX, YY, ZZ);

end