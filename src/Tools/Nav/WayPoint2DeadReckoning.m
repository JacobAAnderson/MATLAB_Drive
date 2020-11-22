% Jacob Anderson
% July 22, 2020

% Convert a waypoint path to a dead reckoning path
% Speed is optional
% Acceptible speed formats:
%     Format            | speed_fmt
%   * meters per second | 'm','mps','m/s'
%   * Knots             | 'knot', 'knots'
%   * Miles per hour    | 'm/h', 'mph'

% Path input must be lon-lat



function [distance, heading, time] = WayPoint2DeadReckoning(path, speed, speed_fmt, dispIO)


% -- Calculate distances form the lat-lon coordinates: s = vdist(lat1,lon1,lat2,lon2) --
for ii = size(path,1)-1 : -1: 1
    
    distance(ii) = vdist(path(ii,2), path(ii,1), path(ii+1, 2), path(ii+1, 1) );
    
    xDists(ii) = vdist(path(ii, 2), path(ii, 1), path(ii, 2),   path(ii+1, 1)) * sign(path(ii+1, 1) - path(ii, 1)); 
    yDists(ii) = vdist(path(ii, 2), path(ii, 1), path(ii+1, 2), path(ii, 1))   * sign(path(ii+1, 2) - path(ii, 2));
    
end

heading = atan2(yDists, xDists);                % Calculate headings

heading = CompassAngle(heading, 'rad', 'deg');  % Conver heading to compass angle


% -- If speed is given, conver to meters per second --
if nargin == 1, return                          % No speed was given, just return the distances
elseif nargin == 2, a = 1;                      % No speed format was given, assume meters per second
else                                            % Conver speed to meters per second 
    switch lower(speed_fmt)
        case {'m','mps','m/s'}, a = 1;          % Format is meter per second
        case {'knot', 'knots'}, a = 1.94384;    % Format is knots
        case {'m/h', 'mph'},    a = 2.23694;    % Format is miles per houre
    end
end


time = distance / (speed*a);                    % Calculate times


% -- Show results --
if nargin == 4 && dispIO
    
    fprintf('\n\n___Waypoint to Dead Reckoning____\nPath:\n')
    disp(path)
    disp(" ")

    disp("Course")

    fprintf('Dist: %f [m], Heading: %f [deg], Duration: %f [s]\n', [distance', heading',time']')
    disp(" ")

end


end


