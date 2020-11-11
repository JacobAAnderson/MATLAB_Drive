

classdef Evolutionalry
    
    properties
        
        Map
        Population
        Start
        MaxDist
        Seeds
    end
    
 
    methods
        
        % Constructor
        function obj = Evolutionalry(map, start, num_Sample, maxDist, seeds)
            
            obj.Map        = map;
            obj.Start      = start;
            obj.MaxDist    = maxDist;
            obj.Population = struct('pop', [], 'area', [], 'info', [], 'score', [], 'dist', [], 'uncert', []);
            obj.Seeds      = seeds;
            
            obj = obj.InitPop(20, num_Sample);
            
        end
        
       
        % Create Path
        function path = InformativePath(obj, epochs, plotIO)
            
            for a = 1: epochs
                
                new_pop = obj.Population;
                old_pop = obj.Population;
            
                for ii = 1: size(obj.Population,2)                    
                    old_path = obj.Population(ii).pop;
                    new_pop(ii)= obj.EvalutatePath( obj.Evolve(old_path) );
                end
                
                old_scores = [obj.Population.score];
                new_scores = [new_pop.score];
                    
                for ii = 1 : size(obj.Population,2) 
                    
                    if min(old_scores) < min(new_scores)
                        ind = find(old_scores == min(old_scores));
                        obj.Population(ii) = old_pop( ind(1));
                        old_scores(ind(1)) = Inf;
                        
                    else
                        ind = find(new_scores == min(new_scores));
                        obj.Population(ii) = new_pop(ind(1));
                        new_scores(ind(1)) = Inf;
                    end
                    
                    
                end
            end
            
            scores = [obj.Population.score];
            path = obj.Population(scores == min(scores)).pop;
            
            if nargin > 2 && plotIO
                obj.Map.PlotMap
                hold on
                plot(path(:,1), path(:,2))
                hold off
            end
            
            
        end
        
    end
    
    
    methods (Access = private)
        
        function obj = InitPop(obj, pop_size, num_samples)
            
            if ~isempty(obj.Seeds)
                s_t = 0;
                for jj = 1: size(obj.Seeds,2), s_t = s_t + obj.Seeds(jj).s; end
            end
            
            for ii = 1: pop_size
                if isempty(obj.Seeds)
                    obj.Population(ii) = obj.EvalutatePath([obj.Start(1:2); obj.Map.RandomPoints(num_samples)]);
                    
                else
                    path = zeros(num_samples,3);
                    start = 1;
                    
                    for jj = 1: size(obj.Seeds,2)
                        
                        mu    = obj.Seeds(jj).mu;
                        sigma = obj.Seeds(jj).sigma;
                        rho   = obj.Seeds(jj).r;
                        s     = obj.Seeds(jj).s;
                        
                        sigma = [sigma(1)^2,   rho * sigma(1) * sigma(2);
                            rho * sigma(1) * sigma(2), sigma(2)^2];
                        
                        num = ceil( s * num_samples / s_t);
                        
                        vals = NaN;
                        while any(isnan(vals))                              % Make Sure the locaitons stay on the map
                            points = mvnrnd(mu,sigma, num);
                            vals = obj.Map.Measure(points);
                        end
                        
                        path( start: start + num-1, 1:2) = points;
                        
                        start = start + num;
                        
                    end
                    
                    obj.Population(ii) = obj.EvalutatePath([obj.Start; path]);
                    
                end
            end
            
            obj.MaxDist = max([obj.Population.dist]);
            
            
        end
        
        
        function path = Evolve(~, path)
            
            ind = [1,1];
            while any(ind ==1), ind = randi(size(path,1), 1,2); end

            a = path(ind(1),:);                                             % Change tour orders
            path(ind(1),:) = path(ind(2),:);
            path(ind(2),:) = a;
            
        end
        
        
        function [info, area, uncert] = Info(obj, path)
            
            ind = zeros(1,3000);
            start = 1;
            for ii = 2: size(path,1)
                idx = obj.FindPointsOnLine(path(ii-1,:), path(ii,:));
                ind(start: start + numel(idx)-1) = idx;
                start = start + numel(idx);
            end
            
            ind(ind ==0) = [];
            
            ind = unique(ind);
            
            ind(isnan(ind)) = [];                                           % Hack !!!!!!!!!!!!
            
            tf = obj.Map.Layers.("Visited")(ind);                           % Eliminate the score from celss that have already been visiteds
            
            if any( isnan(tf))                                              % Check that we stay on the map
                info = Inf;                                                 % Inf value will prevent the path from being selected
                area = Inf;
                uncert = Inf;
                return
            end
            
            ind(logical(tf)) = [];
            
            info = nansum( obj.Map.Map(ind));
            area = obj.Map.GridArea * numel(ind);
            uncert = 0; % sum(obj.Map.Layers.('Uncertainty')(ind));
            
        end
        
        
        function pop = EvalutatePath(obj, path)
            
            % Calculate Path Values
            dist = Distance(path);
            [info, area, uncert] = obj.Info(path);
            
            % Assign values
            pop.pop    = path;
            pop.dist   = dist;
            pop.info   = info;
            pop.area   = area;
            pop.uncert = uncert;
            pop.score  = dist;

            % Multi-objective Optimization
            info = (nansum(obj.Map.Map(:)) - info) / nansum(obj.Map.Map(:)); 
%             area = (obj.Map.Area - area) / obj.Map.Area;
            
%             total_uncert = nansum(obj.Map.Layers.('Uncertainty')(:));
%             uncert = (total_uncert - uncert) / total_uncert;
%             uncert = 0;

%             dist = dist / obj.MaxDist;
%             pop.score  = sqrt(dist*dist + info*info); % + area*area + uncert*uncert);
            
            
        end
        
        
        % Find Points along a line -----------------------------------------------------------------------------------
        function ind = FindPointsOnLine(obj, p1, p2)
            
            m = (p2(2) - p1(2) ) / (p2(1) - p1(1) );    % Line equation
            b = p1(2) - m * p1(1);
            
            x = min(p1(1), p2(1)) : max(p1(1), p2(1));
            
            y = m .* x + b;
            
            ind = obj.Map.state2index([x',y']);
            ind = unique(ind);
            
        end
        
    end
end


function dist = Distance(path)
dist = 0;
for ii = 2: size(path,1)
    diff = path(ii,:) - path(ii-1,:);
    dist = dist + sqrt(diff * diff');
end
end





