


function auv = AUV_Basic( name, wpFile)
% ___ Vehicle kinematics __________________________________________________
% Kinimaticks of a simple car model
% State:  x = [ X; Y; Heading; Speed]
% Control u = [acelleration; string angle]
% Noise   m = [uncertianty in heading; undertainty in speed];

d = 1.5;                                                                    % Wheel base [m]

X = @(x,u,t) [t*x(4)*cos(x(3));                                             % Process Model
              t*x(4)*sin(x(3)); 
              x(4)* tan(u(2))/d; 
                         t*u(1)];

dX = @(x,u,t) [ 1  0  -t*x(4)*sin(x(3))       t*cos(x(3));                  % Proces Model Jacobian
                0  1   t*x(4)*cos(x(3))       t*sin(x(3));
                0  0                  1  x(4)*tan(u(2))/d;
                0  0                  0                 1];    


M = @(x,u,t) [ cos(x(3)), 1;                                                % Process Noise Model
               sin(x(3)), 1; 
                       0, 0; 
                       0, 0];
                          
dM = @(x,u,t) [ 0, 0, -sin(x(3)), 0;
                0, 0,  cos(x(3)), 0;
                0, 0,          0, 0;
                0, 0,          0, 0];
                    
F  = @(x,u,t,m) x + X(x,u,t) + M(x,u,t)*m;                                  % Process Model

rad = @(x) x*pi/180;

% ___ Instanciate the AUV _________________________________________________
auv = R0b0t('AUV', name);                                                   % Instanciate the robot

auv = auv.SetDynamics('f', @(x,u,t,m) F(x,u,t,m) );                         % Add Process Model
auv = auv.SetDynamics('fx_dx', @(x,u,t,m) 1 + dX(x,u,t)+ dM(x,u,t) );       % Add Process Model derivative

auv = auv.SetUncertainty('Speed',  'Normal', -0.05, 0.5);                   % Set uncertainty in state transitions
auv = auv.SetUncertainty('Heading','Normal', rad(5), rad(4));
auv = auv.SetUncertainty('Pitch',  'Normal',  0.0,  5.0);

auv = auv.Add_Sensor("Depth",  'Normal', 0, 0.0597);                        % Add Alitimeter with noise
auv = auv.Add_Sensor("Compass",'Normal', rad(3), rad(6));                   % Add Compass with noise
% auv = auv.Add_Sensor("GPS",    'Normal', 0, 1.5);


% ___ Waypoint Following __________________________________________________
auv.Nav = auv.Nav.Load_Waypoints(wpFile);                                   % Give the ASV a set of waypoints to follow

start = auv.Nav.Wp.Waypoints(1,1:2);                                        % Put the Robot at the First waypoint
auv   = auv.Set_VehicleState(start);

auv.Goal = auv.Nav.Wp.Waypoints(end,1:2);                                   % Put the Last waypoint as the goal  

end

