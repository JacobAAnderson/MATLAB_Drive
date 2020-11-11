% Jacob Anderson
% July 22, 2020

% Convert a waypoint path to a dead reckoning path
% Speed is optional
% Acceptible speed formats:
%     Format            | speed_fmt
%   * meters per second | 'm','mps','m/s'
%   * Knots             | 'knot', 'knots'
%   * Miles per hour    | 'm/h', 'mph'
            



function [distance, heading, time] = WayPoint2DeadReckoning(path, speed, speed_fmt)

fprintf('\n\n___Waypoint to Dead Reckoning____\nPath:\n')
disp(path)
disp(" ")

for ii = size(path,1)-1 : -1: 1
    
    % vdist(lat1,lon1,lat2,lon2)
    
    distance(ii) = vdist(path(ii,2), path(ii,1), path(ii+1, 2), path(ii+1, 1) );
    
    xDists(ii) = vdist(path(ii, 2), path(ii, 1), path(ii, 2),   path(ii+1, 1)) * sign(path(ii+1, 1) - path(ii, 1));
    yDists(ii) = vdist(path(ii, 2), path(ii, 1), path(ii+1, 2), path(ii, 1))   * sign(path(ii+1, 2) - path(ii, 2));
    
end

heading = atan2(yDists, xDists);

heading = CompassAngle(heading, 'rad', 'deg');

disp("Course")
if nargin == 1
    fprintf('Dist: %f [m], Heading: %f [deg]\n', [distance', heading']')
    disp(" ")
    return

elseif nargin == 2, a = 1;
  
else % Conver speed to meters per second
    
    switch lower(speed_fmt)
        
        case {'m','mps','m/s'}, a = 1;
        case {'knot', 'knots'}, a = 1.94384;
        case {'m/h', 'mph'},    a = 2.23694;
            
    end

end

time = distance / (speed*a);


fprintf('Dist: %f [m], Heading: %f [deg], Duration: %f [s]\n', [distance', heading',time']')
disp(" ")

end