

function bathymetry = LoadBathymetry(bathFile), warning on backtrace
bathymetry = load(bathFile);                                                % load Bathymetry map

flds = fields(bathymetry);

if numel(flds) > 1, warning("Bathymetry File Has Too Many Fields"); end

bathymetry = bathymetry.(flds{1});                                          % Get rid of loading structue

end