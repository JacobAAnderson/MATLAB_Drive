%% Side Scan Sonar Shape From Shading

function [zPort,zStar] = BuildHeightMap(intenseImage,  deep, utm, s, maxRange, xOffset, trench)

intenseImage( isnan(intenseImage) ) = 0;

height   = 0.6;

sideLength = size(intenseImage,1)/2-xOffset;


intData(:,:,1) = double( intenseImage( end - sideLength + 1 :end, 1:end ) );    % Port side data
intData(:,:,2) = double( intenseImage( 1:sideLength,         1:end ) );         % Sarboard side data

intData(:,:,2) = flipud(intData(:,:,2));    % Flip the image to mathc the algorithm indexing


rows    = size(intData,1); 
columns = size(intData,2);

Z = zeros(rows, columns,2);
X = zeros(rows, columns,2);

dR = maxRange/(sideLength-1);
ranges = 0: maxRange/(sideLength-1): maxRange;

shadowThresh = 20/255;
maxShadow    = 1;

for sheet = 1:2
    
    for ii = 1:columns
                
        startRow = trench(sheet,ii);
        
        height = ranges(startRow) * cosd(30);
        x      = ranges(startRow) * sind(30);
        xlast  = x;
        dx     = dR * sind(30);
%         Z(1,ii,sheet) = height;
        
        for jj =  startRow : rows-1
            
            x = xlast + dx;
            
            if(intData(jj, ii, sheet ) < shadowThresh)
                continue
            else
                
                psi   = atan2(x, height);                               % Sonar ray angle from horizontal
                theta = acos( sqrt(intData(jj, ii, sheet) ) * 0.7);     % Angle to normal Vector
                alpha = psi - theta;                                    % Inclination of surface plane
                
                
                
                %if(inShadow == 1)
                %    if(lastGood > -1 && (x-xlast) < maxShadow)
                %        Z(j,i) = Z(lastGood,i) - (x-xlast)*tan(lastAlpha);
                %Z(j,i) = Z(lastGood,i) - (x-xlast)/tan(lastAlpha);
                %    end
                %end
                
                if(xlast > -1 && (x-xlast) < maxShadow)
                    Z( jj+1, ii, sheet) = Z( jj, ii, sheet) + (x-xlast)*tan(alpha);   % Calculate next Z value
                    X( jj+1, ii, sheet) = x;
                end
                
                xlast = x;
                x = dR * cos(psi);
                
            end
            
        end
        
    end
    
    y = linspace(0,s,columns);
    
    figure;
    surf( y, X(:,:,sheet), Z(:,:,sheet), 'EdgeColor','none');
    hold on;
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z(m)');
    view(3);
    axis equal
    hold off
    
end

zPort = Z(:,:,1);
zStar = Z(:,:,2);

end



%% Stuff
%% Do 3D Reconstruction
% close all
% clc
% 
% vehicle.X =  0;                         % Vehicle position [metes]
% vehicle.Z = 20;
% 
% 
% range = Data.portSonar(1).range;
% [row, column] = size(portEcho);
% ranges = 0: range/(row-1): range;


% theta = 30;                                     % Sonar angle in degrees up from verticle down [deg]
% Xo = vehicle.X + ranges(ii)   * sind(theta);    % Current echo position
% Zo = vehicle.Z - ranges(ii+1) * cosd(theta);
% 
% 
% % Plot Results
% figure('Name','Triangulation','NumberTitle','off')
% hold on
% plot(vehicle.X, vehicle.Z,'*')
% 
% %scatter3(X(:),Y(:),Z(:),'.')
% viscircles([vehicle.X,vehicle.Z], ranges(ii),'Color','m');
% viscircles([vehicle.X,vehicle.Z], 30,'Color','c');
% axis equal
% 
% % Triangulate next position
% for ii = ii : length(portScanLine) - 1
%     
%     range1 = ranges(ii);                % Sonar range [meters]
%     range2 = ranges(ii+1);
%     
%     Is = portScanLine(ii);              % Sonar and reflectixity charactureistics for Lambert's reflectoion Law
%     % u  = 1;
%     % Ii = 1;
%     % dA = 1;
%     %
%     % p = Is / ( u * Ii * dA);
%     
%     alpha = acosd(2*Is -1) / 2;         % Angle of incidence between the sonar beam and the normal vector of the reflecting surface [deg]
%     alpha2 = acosd(2* portScanLine(ii+1) -1) / 2;
%     beta  = theta - alpha;              % Inclination of surface plan from horizontal [deg]
%     
%     
%     % Calculate triangle betwenn consecutive sonar retruns
%     A = alpha + 90;
%     B = asind(range1 * sind(A) / range2 );
%     C = 180 - ( A + B );
%     
%     % Third side of the triangle
%     c = sqrt(range1*range1 + range2*range2 - 2 * range1 * range2 * cosd(C));
%     
%     
%     deltaX = c * cosd(beta);            % Change in X and Z
%     deltaZ = c * sind(beta);
%     
%     X(ii ) = Xo + deltaX;               % Next echo point
%     Z(ii ) = Zo + deltaZ;
%     
%     fprintf('\n %f, %f,  %f,  %f', alpha, beta, C, theta)
%     
%     Xo = X( ii );
%     Zo = Z( ii );
%     theta = theta + C;
%     
%     plot(X(ii) ,Z(ii),'.','Color','g')
%     
%     
% end


% plot(X,Z,'.','Color','g')
% plot(X,Z,'*','Color','r')

% axis equal
% hold off
% drawnow

%end


% Ii = 1:1:length(portScanLine);
% alpha = acosd(2 * portScanLine -1) / 2; 
% 
% figure('Name','Alpha')
% plot(Ii, alpha)

