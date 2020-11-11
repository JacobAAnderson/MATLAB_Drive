%% Get Dearreckoning and SLAM path Structure
%  Jacob Anderson
%  8/1/2019

function path = Get_drPath_struct(emData, showPath, Geotiff)

[path.dr(:,1), path.dr(:,2)] = Vel_head2Path( emData.speed, emData.heading, emData.timestamp);
path.dr(:,1) = path.dr(:,1) + emData.utm(1,1);
path.dr(:,2) = path.dr(:,2) + emData.utm(1,2);

% Get the differance between the dead reckoning path and the gps path to adjust the bathymetry readings
path.diff(:,1) = path.dr(:,1) - emData.utm(:,1);
path.diff(:,2) = path.dr(:,2) - emData.utm(:,2);

% What to process ----------------------------------------------------
path.corrected(:,1) = path.dr(:,1);
path.corrected(:,2) = path.dr(:,2);

% Location and depth of altimeter readings
path.altimeter(:,1) =  emData.altimiter(:, 1) + path.diff(:,1);
path.altimeter(:,2) =  emData.altimiter(:, 2) + path.diff(:,2);

% Deturmin uncertainties ---------------------------------------------------------------------
[path.err(:,1), path.err(:,2)] = PathError(emData.utm(:,1), emData.utm(:,2), emData.speed, emData.heading, emData.timestamp);



if nargin > 2 && showPath
    PlotPathes(emData, path, Geotiff)
    
elseif nargin > 1 && showPath
    PlotPathes(emData, path)
    
else
    return
end

end


%% Calculate Speed and Heading Error
function[errTheta, errSpeed] = GetVelocityErrors(gpsX, gpsY, speed, heading, timestamp)

% Get Path from Speed and heading
timeStep = seconds(timestamp(2:end) - timestamp(1: end-1));

deltaX_dr = [0; speed(1:end-1) .* timeStep .* sin( heading(1:end-1) )];
deltaY_dr = [0; speed(1:end-1) .* timeStep .* cos( heading(1:end-1) )];

deltaX_gps = [0; gpsX(2:end) - gpsX(1:end-1)];
deltaY_gps = [0; gpsY(2:end) - gpsY(1:end-1)];

% Error between deadreckoning and GPS
errX = deltaX_dr - deltaX_gps;
errY = deltaY_dr - deltaY_gps;


% Deturmin uncertainties ---------------------------------------------------------------------
sigtrig = sin(heading).^2 - cos(heading).^2;

sigX_dX = errX ./ deltaX_dr;
sigY_dY = errY ./ deltaY_dr;

errTheta = abs(heading) .* sqrt( abs((sigY_dY .* sigY_dY  - sigX_dX .* sigX_dX) ./ sigtrig ));
errSpeed = abs(speed)   .* sqrt( abs((tan(heading ).^2 .*  sigX_dX .^2 - sigY_dY.^2) ./ ( tan(heading ).^2 - 1)));

errTheta(isinf(errTheta)) = [];
errSpeed(isinf(errSpeed)) = [];

errTheta = nanmean(errTheta);
errSpeed = nanmean(errSpeed);

fprintf("\n\nTheta Err: %f\tSpeed Err: %f\n\n",errTheta, errSpeed)

end


%% Path From Deadreckoning
function [pathX, pathY] = Vel_head2Path( speed, heading, timestamp)

timeStep = seconds(timestamp(2:end) - timestamp(1: end-1));

deltaX_dr = speed(1:end-1) .* timeStep .* sin( heading(1:end-1) );
deltaY_dr = speed(1:end-1) .* timeStep .* cos( heading(1:end-1) );

pathX = [0; cumsum(deltaX_dr)];
pathY = [0; cumsum(deltaY_dr)];

end


%% Path Error
function [errX, errY] = PathError(gpsX, gpsY, speed, heading, timestamp)

[pathX, pathY] = Vel_head2Path( speed, heading, timestamp);

% [errTheta, errSpeed] = GetVelocityErrors(gpsX, gpsY, speed, heading, timestamp);
errTheta = 0.08;    % Heading Error in rads 
errSpeed = 0.25;    % Speed Error in m/s

deltaX_dr = [0; pathX(2:end) - pathX(1:end-1)];
deltaY_dr = [0; pathY(2:end) - pathY(1:end-1)];

sigV_v = errSpeed ./ speed;
sigT_t = errTheta ./ heading;

sigV_v(isinf(sigV_v)) = 0;
sigT_t(isinf(sigT_t)) = 0;

errX = abs(deltaX_dr) .* sqrt( sigV_v.^2 + ( cos(heading) .* sigT_t ).^2 );
errY = abs(deltaY_dr) .* sqrt( sigV_v.^2 + ( sin(heading) .* sigT_t ).^2 );

errX = cumsum(errX);
errY = cumsum(errY);

end


%% Plot Processed Data
function PlotPathes(emData, path, Geotiff)

figure('Name','Pathes', 'NumberTitle','off')
hold on

if nargin > 2
    geoshow(Geotiff.image, Geotiff.data)
end
    
p1 = plot(emData.utm(:,1),     emData.utm(:,2),    'b');
p2 = plot(path.dr(:,1),        path.dr(:,2),       'r');
p3 = plot(path.corrected(:,1), path.corrected(:,2),'g');

scatter(path.dr(1,1), path.dr(1,2),'d','Cdata', [0,0,0])

scatter(path.dr(end,1), path.dr(end,2),'*','Cdata', [1,0,0])
scatter(path.dr(end,1), path.dr(end,2),'o','Cdata', [0,0,0])

scatter(emData.utm(end,1), emData.utm(end,2),'*','Cdata', [0,0,1])
scatter(emData.utm(end,1), emData.utm(end,2),'o','Cdata', [0,0,0])

scatter(path.corrected(end,1), path.corrected(end,2),'*','Cdata', [0,1,0])
scatter(path.corrected(end,1), path.corrected(end,2),'o','Cdata', [0,0,0])

hold off

legend([p1,p2,p3],{'GPS','Dead Reckoning','Bathymetry & KF'})
xtickformat('%g')
ytickformat('%g')

end
