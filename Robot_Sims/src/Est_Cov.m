% Dead Reckoning Calculation
% Jacob Anderson
% 2/11/2020

function Cov = Est_Cov(speed, sig_s, heading, sig_h, dt, heading_fmt, coor_fram)

if nargin > 5 && strcmpi(heading_fmt, 'deg'), heading = heading * pi/180; end

if speed == 0 
    Cov = [0.1, 0; 0, 0.1];
    return
end

sig_s2 = (sig_s./speed).*(sig_s./speed);

sig_x = abs(dt) .* sqrt( sig_s2 + (sin(heading).* sig_h./heading) .* (sin(heading).* sig_h./heading) );
sig_y = abs(dt) .* sqrt( sig_s2 + (cos(heading).* sig_h./heading) .* (cos(heading).* sig_h./heading) );

if nargin > 6 && strcmpi(coor_fram, 'compass')
    
    sig_xy = -0.8 * cos( 2*heading ) .* sig_x .* sig_y;
    Cov = [sig_y*sig_y, sig_xy; sig_xy, sig_x*sig_x];
    
else
    sig_xy = -0.8 * sin( 2*heading ) .* sig_x .* sig_y;
    Cov = [sig_x*sig_x, sig_xy; sig_xy, sig_y*sig_y];
    
end

end
