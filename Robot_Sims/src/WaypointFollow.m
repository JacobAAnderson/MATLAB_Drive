% Waypoint Follower
% Jacob Anderson
% 2/12/2020

% Returns normalized vehicle speed --> Multiply by vehicle's max speed
% heading in radians from positive X in x-y-z(up) frames

function [speed, bearing, wp_idx, done] = WaypointFollow(X, path, r, wp_idx)

done = false;

goal = path(wp_idx,:);

dist = goal(1:2) - X(1:2)';                                                 % Distance between robot and goal in cartesian spcae

speed = tanh(vecnorm(dist));                                                % Normalized Speed

bearing = atan2(dist(2), dist(1));                                          % bearing to the goal


if sqrt(dist' * dist) < r, wp_idx = wp_idx +1; end                          % Check for waypoint achival
if wp_idx > size(path,1), done = true; wp_idx = size(path,1);  end          % Check if the path is completed

end