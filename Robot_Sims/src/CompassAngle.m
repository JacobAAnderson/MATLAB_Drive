% Convert to Compass angle
% Jacob Anderson
% 2/11/2020


function heading = CompassAngle(heading, in_fmt, out_fmt)

if strcmpi(in_fmt,'rad'), heading = heading*180/pi;                         % Conver heading to degrees
elseif strcmpi(in_fmt,'deg')
else, warning('Un-reckonized input format, no converstions have been applied')
end

if heading < 90, heading = 90 - heading;                                    % Conver to Compass angle
else, heading = 450 - heading;
end

if strcmpi(out_fmt,'rad'), heading = heading*pi/180;                        % Conver heading to degrees
elseif strcmpi(out_fmt,'deg')
else, warning('Un-reckonized output format, no converstions have been applied')
end

end