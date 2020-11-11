

classdef WaterFeature
    
    properties
        Map
        Grad
        Features
        Model
    end
    
    properties(Access = private)
        Gaussians
        XX
        YY
        indx
        noise
        spreading
    end
    
    
    methods
        
        function obj = WaterFeature(XX, YY, map, name)
            
            if nargin > 3, obj.Map = Mapp(XX, YY, map, name);
            else, obj.Map = Mapp(XX, YY, map);
            end
            
            obj.Map.Map(~isnan(map)) = 0;             % Set Map to 0
            
            % Create local referance frame ------------------
            [r,c] = size(map);
            
            x = linspace(-10, 10, c);
            y = linspace(-10, 10, r);
            
            [obj.XX ,obj.YY] = meshgrid(x,y);
            
            % Gaussians to model the water features --------
            obj.Gaussians = struct('Mu', [], 'Sigma',[], 'Rho', [], 'Scale', []);
            obj.indx = 1;
            
            % Set noise levels to 0
            obj.noise.Compass   = makedist('Normal','mu',0,'sigma',0);      
            obj.noise.Speed     = makedist('Normal','mu',0,'sigma',0);      
            obj.noise.Spreading = makedist('Normal','mu',0,'sigma',0);
        end
        
        
        % ____ Set Up _____________________________________________________
        % Set vehicle transition model
        function obj = SetUncertainty(obj, type, dist, mu, sig)
            obj.noise.(type) = makedist(dist,'mu',mu,'sigma',sig);
        end
     
        
        % Add pices to the feature
        function obj = Add_Gaussian(obj, mu, sigma, rho, s)
            
            sigma = [sigma(1)^2,   rho * sigma(1) * sigma(2);
                rho * sigma(1) * sigma(2), sigma(2)^2];
            
            obj.Gaussians(obj.indx).Mu = mu;
            obj.Gaussians(obj.indx).Sigma = sigma;
            obj.Gaussians(obj.indx).Scale = s;
            obj.indx = obj.indx + 1;
            
        end

        
        % ____ Processes __________________________________________________
        % Move the Feature
        function obj = Move(obj, speed, theta, dt)
            
            theta = theta + random(obj.noise.Compass, 1);
            theta(theta > 360) = theta(theta > 360) - 360;
            theta(theta < 0)   = theta(theta < 0)   + 360;
            
            speed = speed + random(obj.noise.Speed,1);
            
            dx = speed * dt * cosd(theta);
            dy = speed * dt * sind(theta);
            
            for ii = 1: obj.indx - 1
                
                sig = random(obj.noise.Spreading,3);
                mu    = obj.Gaussians(ii).Mu;
                obj.Gaussians(ii).Mu = [mu(1) + dx + sig(1), mu(2)+dy+ sig(2)];
            end
        end
        
        
        % ____ Data Processing ____________________________________________
        function obj = MakeFromData(obj, X, y)
            
            obj.XX = obj.Map.Easting;
            obj.YY = obj.Map.Northing;
            
            obj.Model.GP = fitrgp(X,y);
            
            Xnew = [obj.Map.Easting(:),  obj.Map.Northing(:)];
            sz = size(obj.Map.Map);
            
            [ypred,ysd] = predict(obj.Model.GP,Xnew);
            
            ypred = reshape(ypred, sz);
            ysd   = reshape(ysd,   sz);
            
            ypred(isnan(obj.Map.Map)) = NaN;
            ysd(isnan(obj.Map.Map)) = NaN;
            
            obj.Map.Map = ypred;                                                % GP Model Becomes the new Map
            obj.Map = obj.Map.Add_Layer('Uncertainty', ysd);
            
            obj.Model.ypred = ypred;
            obj.Model.ysd = ysd;
            
        end
        
        
        function obj = Make3DModelFromData(obj, X, y, bathy)
            
            obj.Model.GP = fitrgp(X,y);
            
            x = linspace(min(obj.Map.Easting(:)),  max(obj.Map.Easting(:)),  size(obj.Map.Easting,2));
            y = linspace(min(obj.Map.Northing(:)), max(obj.Map.Northing(:)), size(obj.Map.Easting,1));
            z = 0: -0.25: floor(min(-bathy(:)));
            
            [xx, yy, zz] = meshgrid(x, y, z);
            
            for ii = 1: size(xx,3)
                
               x = xx(:,:,ii);
               y = yy(:,:,ii);
               z = zz(:,:,ii);
               
               tf = z < -bathy;
               
               x(isnan(bathy) | tf) = NaN;
               y(isnan(bathy) | tf) = NaN;
               z(isnan(bathy) | tf) = NaN;
               
               xx(:,:,ii) = x;
               yy(:,:,ii) = y;
               zz(:,:,ii) = z;
               
            end
            
            Xnew = [xx(~isnan(xx)), yy(~isnan(xx)), zz(~isnan(xx))];
            
            [ypred,ysd] = predict(obj.Model.GP,Xnew);
            
            obj.Model.xx = xx(~isnan(xx));
            obj.Model.yy = yy(~isnan(xx));
            obj.Model.zz = zz(~isnan(xx));
            obj.Model.ypred = ypred;
            obj.Model.ysd = ysd;
            
        end
 
        
        % Build Feature Model
        function obj = MakeMap(obj)
            
            map = zeros(size(obj.Map.Map));
            
            for ii = 1: obj.indx - 1
                
                mu    = obj.Gaussians(ii).Mu;
                sigma = obj.Gaussians(ii).Sigma;
                s     = obj.Gaussians(ii).Scale;
                
                Y = mvnpdf([obj.XX(:), obj.YY(:)], mu, sigma);
                
                map = map + s* reshape(Y, size(obj.XX));
                
            end
            
            map(isnan(obj.Map.Map)) = NaN;
            
            obj.Map.Map = map;
            
        end
        
        
        % Build Feature Model
        function obj = MakeMap_world(obj)
            
            map = zeros(size(obj.Map.Map));
            
            for ii = 1: obj.indx - 1
                
                mu    = obj.Gaussians(ii).Mu;
                sigma = obj.Gaussians(ii).Sigma;
                s     = obj.Gaussians(ii).Scale;
                
                Y = mvnpdf([obj.Map.Easting(:), obj.Map.Northing(:)], mu, sigma);
                
                map = map + s* reshape(Y, size(obj.XX));
                
            end
            
            map(isnan(obj.Map.Map)) = NaN;
            
            obj.Map.Map = map;
            
        end
  
        
        % Find grdainat of the feature and the max/mins of the second derivative
        function obj = Find_Gradiant(obj)
            
            % Calculate the second derivative of the scalar field
            [FX,FY] = gradient(obj.Map.Map);
            
            FXX = gradient(FX);
            [~,FYY] = gradient(FY);
            
            grad =  (FXX + FYY);                    % Second deriviative of the scaler field
            
            % Find local Maxima and minima
            TFc = islocalmin(grad,1);
            TFr = islocalmin(grad,2);
            
            TFmin = TFc & TFr;                      % Indicies of local minima
            
            TFc = islocalmax(grad,1);
            TFr = islocalmax(grad,2);
            
            TFmax = TFc & TFr;                      % Indicies of local maxima
            
            % Assign Values to model
            obj.Grad.grad = grad;
            obj.Grad.max_points = [obj.XX(TFmax), obj.YY(TFmax), obj.Map.Map(TFmax)];
            obj.Grad.min_points = [obj.XX(TFmin), obj.YY(TFmin), obj.Map.Map(TFmin)];
            
            
            % Cluster points by height to eliminate outliers
            idx = kmeans(obj.Grad.min_points(:,3), 2);
            
            c1 = mean(obj.Grad.min_points(idx==1, 3));
            c2 = mean(obj.Grad.min_points(idx==2, 3));
            
            if c1 > c2, obj.Grad.min_points = obj.Grad.min_points(idx==1, :);
            else, obj.Grad.min_points = obj.Grad.min_points(idx==2, :);
            end
            
        end
        
        
        % Find Standard deviation of the features
        function obj = Fit_Model(obj)
            
            % Get Sigma Estimates
            ptCloud = pointCloud(obj.Grad.min_points);
            testPoints = obj.Grad.max_points;
            
            match   = zeros(size(testPoints));
            
            for ii = 1: size(testPoints,1)                                     % Conect maxima to minima
                [indices,~] = findNearestNeighbors(ptCloud, testPoints(ii,:), 1);
                match(ii,:) = obj.Grad.min_points(indices,:);
            end
            
            dists = abs(testPoints - match);
            
            [obj.Grad.min_points, ~, n] = unique(match, 'rows');
            
            A = cat(1, obj.Grad.min_points(:,3));
            
            % Initialize Gausians
            for ii = 1: size(obj.Grad.min_points,1)
                
                sigma = mean( dists(n ==ii,:), 1)/2;
                sigma = sigma(1:2);
                
                if all(sigma == 0)
                    warning("Its ALLL WRONG !!!")
                elseif any(sigma == 0)                                        % Keep Sigma Positive defined
                    sigma(sigma == 0) = sigma(sigma ~= 0) / 2; 
                end
                
                obj.Features.property(ii).mu    = obj.Grad.min_points(ii,1:2);
                obj.Features.property(ii).sigma = sigma(1:2);
                obj.Features.property(ii).s     = obj.Grad.min_points(ii,3);
                obj.Features.property(ii).r     = 0;
                
            end
            
            obj = obj.Make_Estimate;
            %             fig = obj.PlotModel;
            
            % Fit the Data
            for gaus = 1 : size(obj.Grad.min_points,1)
                
                for rho = -0.8: 0.1: 0.8
                    obj.Features.property(gaus).r = rho;
                    
                    % Deblend Gaussians -----------------------------------
                    B = eye(size(obj.Grad.min_points,1));
                    
                    for ii = 1: size(B,1)
                        for jj = 1: size(B,2)
                            
                            if ii == jj, continue, end
                            
                            sig = obj.Features.property(jj).sigma;
                            mu  = abs(obj.Features.property(ii).mu - obj.Features.property(jj).mu);
                            r   = obj.Features.property(jj).r;
                            
                            u = 2 * r * (mu * mu') / (sig * sig');
                            mu = (mu .* mu)./sig;
                            
                            B(ii,jj) = exp(-1/(2* (1 - r*r)) * ( mu*mu' + u));
                            
                        end
                    end
                    
                    Ar = B\A;
                    
                    % Adjust Scales
                    for ii = 1: size(obj.Grad.min_points,1)
                        
                        obj.Features.property(ii).s = Ar(ii);
                        
                    end
                    
                    obj = obj.Make_Estimate;
                    
                    %                     fig = obj.PlotModel([0,90], fig);
                    
                    D = abs(obj.Map.Map - obj.Features.Model).^2;
                    mse = nansum(D(:))/ sum(~isnan(obj.Map.Map(:)));
%                     mse = MSError(D, obj.Map.Map);
                    
                    if rho == -0.8
                        min_MSE = mse;
                        model = obj.Features;
                        
                    elseif mse < min_MSE
                        
                        min_MSE = mse;
                        model = obj.Features;
                        
                    end
                    
                end
                
                obj.Features = model;
                
            end
            
%             fprintf("\n\nMin MSE: %f\n\n", min_MSE)
            obj.Features = model;
            %             obj.PlotModel
            
        end
        
        
        % Genrate Estimate from the model
        function obj = Make_Estimate(obj)
            
            obj.Features.Model = zeros(size(obj.Map.Map));
            
            for ii = 1: size(obj.Features.property,2)
                
                mu  = obj.Features.property(ii).mu(1:2);
                sig = obj.Features.property(ii).sigma;
                s   = obj.Features.property(ii).s;
                r   = obj.Features.property(ii).r;
                
                sigma = [sig(1)*sig(1), r*sig(1)*sig(2);
                    r*sig(1)*sig(2), sig(2)*sig(2)];
                
                try
                    Y = mvnpdf([obj.XX(:), obj.YY(:)], mu, sigma);
                catch E
                    warning(E.message)
                    disp(mu)
                    disp("-----------------")
                    disp(sig)
                    disp("-----------------")
                    disp(r)
                    continue
                end
                
                Y = Y / nanmax(Y(:));
                
                obj.Features.Model =  obj.Features.Model + s* reshape(Y, size(obj.XX));
                
            end
            
            obj.Features.Model(isnan(obj.Map.Map)) = NaN;
            
        end
        
        
        
        % ____ Potting ____________________________________________________
        % Plot the feature model
        function fig = PlotMap(obj, v, fig)
            
            if nargin > 2, figure(fig),
            else, fig = figure('Name', 'Water Featue', 'numbertitle', 'off');
            end
            
            levels = linspace(nanmin(obj.Map.Map(:)), nanmax(obj.Map.Map(:)), 10);
            
            surf(obj.XX, obj.YY, obj.Map.Map, 'EdgeColor', 'none' )
            colorbar
            hold on
            contour3(obj.XX, obj.YY, obj.Map.Map, levels, 'k')
            
            hold off
            xlabel('X')
            ylabel('Y')
            
            if nargin > 1, view(v(1), v(2))
            else, view(0, 90)
            end
            drawnow
        end
        
        
        % Plot the Gradiant of 2nd derivative
        function fig = PlotGrad(obj, v, fig)
            
            if nargin > 2, figure(fig), hold on;
            else, fig = figure('Name', 'Grad', 'numbertitle', 'off');
            end
            
            levels = linspace(nanmin(obj.Grad.grad(:)), nanmax(obj.Grad.grad(:)), 10);
            
            surf(obj.XX, obj.YY, obj.Grad.grad, 'EdgeColor', 'none' )
            colorbar
            hold on
            contour3(obj.XX, obj.YY, obj.Grad.grad, levels, 'k')
            hold off
            xlabel('X')
            ylabel('Y')
            
            if nargin > 1, view(v(1), v(2))
            else, view(0, 90)
            end
            drawnow
            
        end
        
        
        % Plot Mak min of the gradiant
        function fig = Plot_MaxMin(obj, v, fig)
            
            if nargin > 2, figure(fig),
            else, fig = figure('Name', 'Max-Min', 'numbertitle', 'off');
            end
            
            maxPoints = obj.Grad.max_points;
            minPoints = obj.Grad.min_points;
            
            hold on
            
            s1 = scatter3(minPoints(:,1), minPoints(:,2), minPoints(:,3), '*b');
            scatter3(minPoints(:,1), minPoints(:,2), minPoints(:,3), 'ob')
            
            s2 = scatter3(maxPoints(:,1), maxPoints(:,2), maxPoints(:,3), '*r');
            scatter3(maxPoints(:,1), maxPoints(:,2), maxPoints(:,3), 'or')
            
            hold off
            
            legend([s1,s2],{"Min", "Max"})
            
            if nargin > 1, view(v(1), v(2))
            else, view(0, 90)
            end
            drawnow
            
        end
        
        
        % Plot the feature model
        function fig = PlotModel(obj, v, fig)
            
            if nargin > 2, figure(fig),
            else, fig = figure('Name', 'Water Featue Estimate', 'numbertitle', 'off');
            end
            
            levels = linspace(nanmin(obj.Features.Model(:)), nanmax(obj.Features.Model(:)), 10);
            
            surf(obj.XX, obj.YY, obj.Features.Model, 'EdgeColor', 'none' )
            colorbar
            hold on
            contour3(obj.XX, obj.YY, obj.Features.Model, levels, 'k')
            hold off
            xlabel('X')
            ylabel('Y')
            
            if nargin > 1, view(v(1), v(2))
            else, view(0, 90)
            end
            drawnow
        end
        
        
        function fig = PlotGP_Prediction(obj,fig)
            
            im = 10 * obj.Model.ysd/nanmax(obj.Model.ysd(:));
            
            im(isnan(im)) = 0;
            
            im = round(im);
            
            map = parula(10);
            rgbImage = ind2rgb(im, map);
            
            levels = linspace(nanmin(obj.Model.ypred(:)), nanmax(obj.Model.ypred(:)), 10);

            
            if nargin > 1, figure(fig), hold on,
            else, fig = figure('Name', 'Gaussian Process Prediction', 'numbertitle', 'off');
            end
            
            surf( obj.Map.Easting,  obj.Map.Northing, obj.Model.ypred, rgbImage);
            shading interp
            h = colorbar;
            title(h,'Standard Deviation')
            hold on
            
            contour3(obj.Map.Easting,  obj.Map.Northing, obj.Model.ypred, levels, 'r');
            
            
            view(23,57)
            
            xlabel('X')
            ylabel('Y')
            hold off
            
        end
        
        
        function Plot3DModel(obj, bathy)
            
            
            alphas = 0 : 0.01: 1;
            
            [Y,~] = discretize(obj.Model.ypred,numel(alphas)); 
            
            figure('name','3D Model')
            surf(bathy.Easting, bathy.Northing, -bathy.Map, 'EdgeColor', 'none' )
            hold on
            
            for ii = 1:numel(alphas)
                s = scatter3(obj.Model.xx(Y ==ii) ,obj.Model.yy(Y ==ii) , obj.Model.zz(Y ==ii) ,'.');
                s.CData = [0.0588    1.0000    1.0000];
                s.MarkerEdgeAlpha = alphas(ii);
                s.MarkerFaceAlpha = alphas(ii);
            end
            
        end
    end
      
end
