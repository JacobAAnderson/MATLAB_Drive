% Dead Reckoning Calculation
% Jacob Anderson
% 2/11/2020

function x = DeadReckoning(x, speed, heading, dt, heading_fmt, coor_fram)

if nargin > 4 && strcmpi(heading_fmt, 'deg'), heading = heading * pi/180; end

if nargin > 5 && strcmpi(coor_fram, 'compass')
    x = x + speed * dt .* [sin(heading), cos(heading)];
else  
    x = x + speed * dt .* [cos(heading), sin(heading)];
end

end
