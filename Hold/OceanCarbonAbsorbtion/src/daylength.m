function dl=daylength(yd,lat)
%
% calculates daylength.  required inputs are day of year (yd)
% and latitude (lat)
%
%
dec=23.5.*cos((2.*pi.*(yd - 172))./365);	%%%declination of sun
dl=acosd(-tand(lat).*tand(dec))./15.*2;		%%%sunrise hour angle / 15deg/hr x 2 (sr +ss)
dl=real(dl);
return
