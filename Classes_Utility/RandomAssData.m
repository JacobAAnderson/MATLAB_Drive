

classdef RandomAssData
    
    properties
        Data
        Stats
        Idx
        
    end
    
    methods
        
        % Create Data Class
        function obj = RandomAssData(sz, types, params), warning on backtrace
            
            if nargin == 1
                types = 'Base';
                params = 'Err';
            elseif nargin == 2
                params = 'Err';
            end
            
            for t = types
                for p = params
                    
                    if strcmpi(p{1}, 'time')
                        obj.Data.(t{1}).(p{1}) = NaT(sz);
                    else
                        obj.Data.(t{1}).(p{1}) = NaN(sz);
                    end
                    
                end
                obj.Idx.(t{1}) = 1;
            end
            
        end
        
        
        % Add Data Points
        function obj = Add_Result(obj, iter, type, varargin)
            
                
            if size(varargin,2) == size(fieldnames(obj.Data.(type)),1)
                
                for ii = [ fieldnames(obj.Data.(type))'; varargin ]
                    
                    obj.Data.(type).(ii{1})( obj.Idx.(type), iter ) = ii{2};
                end
                
                obj.Idx.(type) = obj.Idx.(type) +1;
                
            else
                warning("Number of Input Arguments Does Not Match Number of Data Parameters")
            end
        end
        
        
        
        function obj = AddRad(obj, rAd)
            
            flds = fields(rAd.Data);
            
            for ii = 1:numel(flds)
                
                if isfield(obj.Data, flds{ii})
                    
                    types = fields(rAd.Data.(flds{ii}));
                    
                    for jj = 1: numel(types)
                        
                        if isfield(obj.Data.(flds{ii}), types{jj})
                            start = sum(~isnan( obj.Data.(flds{ii}).(types{jj})(1,:)) ) + 1;
                            
                            step = size(rAd.Data.(flds{ii}).(types{jj}), 2) - 1;
                            
                            obj.Data.(flds{ii}).(types{jj})(:,start: start + step) = rAd.Data.(flds{ii}).(types{jj});
                        else
                            obj.Data.(flds{ii}).(types{jj}) = rAd.Data.(flds{ii}).(types{jj});
                        end
                    end
                else
                    obj.Data.(flds{ii}) = rAd.Data.(flds{ii});
                end
            end
            
        end
        
        
        
        % Reset All Varialb eindecies
        function obj = ResetIndex(obj)
            
            for f = fields(obj.Idx)',  obj.Idx.(f{1}) = 1;  end
            
        end
        
        
        % Do Statisitce on Data
        function obj = Eval_Stats(obj, x_param, dispIO)
            
            cal_area = true;
            
            if nargin < 2 || isempty(x_param), cal_area = false; end
            if nargin < 3, dispIO = false; end
            
            if dispIO, fprintf("_________________________________________\nData Stats\n\n"); end
            
            
            for t_ = fieldnames(obj.Data)', t = t_{1};                      % Iterate through types and take t out of is cell
                
                if dispIO, fprintf(" ---- %s ----\n", t), end
                
                for p_ = fieldnames(obj.Data.(t))', p = p_{1};              % Iterate throught parameters and take p out of its cell
                    
                    if isdatetime(obj.Data.(t).(p)(1,:))
                        n = sum( ~isnat( obj.Data.(t).(p)(:,1) ) );         % Get sample size to calculate standared error of the mean
                        
                        if n == 0, continue, end                                % If data has not been added to this parameter, then skip it
                        
                        % Find Max Error
                        d_t = datenum(obj.Data.(t).(p));
                        
                        obj.Stats.(t).(p).Max = datetime( nanmax( d_t, [], 2),          'ConvertFrom', 'datenum' );
                        obj.Stats.(t).(p).Ave = datetime( nanmean(d_t,     2),          'ConvertFrom', 'datenum' );
                        obj.Stats.(t).(p).STD = datetime( nanstd( d_t, [], 2),          'ConvertFrom', 'datenum' );
                        obj.Stats.(t).(p).SEM = datetime( nanstd( d_t, [], 2)./sqrt(n), 'ConvertFrom', 'datenum' );
                        obj.Stats.(t).(p).Min = datetime( nanmin( d_t, [], 2),          'ConvertFrom', 'datenum' );
                        
                        
                    else
                        n = sum( ~isnan( obj.Data.(t).(p)(:,1) ) );
                        
                        if n == 0, continue, end                                % If data has not been added to this parameter, then skip it
                        
                        % Find Max Error
                        obj.Stats.(t).(p).Max = nanmax( obj.Data.(t).(p), [], 2);
                        obj.Stats.(t).(p).Ave = nanmean(obj.Data.(t).(p),     2);
                        obj.Stats.(t).(p).STD = nanstd( obj.Data.(t).(p), [], 2);
                        obj.Stats.(t).(p).SEM = nanstd( obj.Data.(t).(p), [], 2)./sqrt(n);
                        obj.Stats.(t).(p).Min = nanmin( obj.Data.(t).(p), [], 2);
                        
                        
                    end
                    %                     % Find Average Run Length
                    %                     in = ~isnan( obj.Data.(t).(p) );
                    %                     in = sum(in,1);
                    %                     s  = std(in);
                    %                     in = mean(in);
                    %                     in = round(in-s);
                    %
                    %                     % Truncate states to the average length
                    %                     obj.Stats.(t).(p).Max(in : end) = [];
                    %                     obj.Stats.(t).(p).Ave(in : end) = [];
                    %                     obj.Stats.(t).(p).STD(in : end) = [];
                    %                     obj.Stats.(t).(p).SEM(in : end) = [];
                    %                     obj.Stats.(t).(p).Min(in : end) = [];
                    
                    if dispIO
                        fprintf("\t%s --> Max: %f\n",        p, obj.Stats.(t).(p).Max)
                        fprintf("\t%s --> Mav: %f %s  %f\n", p, obj.Stats.(t).(p).Ave, char(177), obj.Stats.(t).(p).STD )
                        fprintf("\t%s --> Min: %f\n\n",      p, obj.Stats.(t).(p).Min)
                    end
                end
                
            end
            
           
            if ~ cal_area, return, end 
            
            
            % ---- Get Areas ---
            for t_ = fieldnames(obj.Data)', t = t_{1};                      % Iterate through types and take t out of is cell
                
                for p_ = fieldnames(obj.Data.(t))', p = p_{1};              % Iterate throught parameters and take p out of its cell
                    
                    if strcmp(p, x_param), continue, end                    % Don't calculate area with respect to yourself
                    
                    xData = obj.Stats.(t).(x_param).Ave;
                    yData = obj.Stats.(t).(p).Ave;                          
                    
                    
                    if ~ ( size(xData,2) == 1 && size(yData,2) == 1 )
                        warning('x and y data have differnt sizes')
                    end
                    
                    
                    diffs = xData(2:end) - xData(1:end-1);
                    
                    
                    if isduration(diffs)
                        diffs = seconds(diffs);
                    end
                    
                    
                    area = nansum( [0; diffs] .* yData );
                    
                    obj.Stats.(t).(p).Area_Ave = area;
                end
                
            end
            
            
            
        end
        
        
        % Write Stats to a table
        function Write_Stats(obj, filename, param, fm)
            
            [path, name, ect] = fileparts(filename);
            
            switch ect
                
                case '.txt'
                    pm = char(177);
                    
                case '.tex'
                    pm = '$\pm$';
                otherwise
                    warning("Unrecognized File Type %s", ect)
                    return
            end
            
            
            types = fieldnames(obj.Data);
            env = cell(3, size(types,1));
            
            if nargin < 4, fm = '%3.2f'; end
            
            for ii = 1: size(types,1)
                
                env{1, ii} = sprintf(fm, obj.Stats.(types{ii}).(param).Max);
                env{2, ii} = sprintf([fm, ' %s ' fm], obj.Stats.(types{ii}).(param).Ave, pm, obj.Stats.(types{ii}).(param).STD);
                env{3, ii} = sprintf(fm, obj.Stats.(types{ii}).(param).Min);
                
            end
            
            stats = {'Max'; 'Ave'; 'Min'};
            
            T = array2table(env, 'RowNames', stats, 'VariableNames', upper(types));
            
            filename = fullfile(path, [name,'_',param, ect]);
            
            switch ect
                
                case '.txt'
                    
                case '.tex'
                    table2latex(T, filename)
                otherwise
                    warning("Unrecognized File Type %s", ect)
            end
            
            
        end
        
        
        % Plot the data
        function PlotData(obj, type, param, vs)
            
            figure('name', param, 'numbertitle', 'off')
            
            n = numel(type);
            
            p(n) = gobjects;
            
            for ii = 1: n
                x_data = obj.Data.(type{ii}).(vs);
                y_data = obj.Data.(type{ii}).(param);
                
                p(ii) = plot(x_data, y_data);
                hold on
            end
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting

            legend(p, type, 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 16, 'Box','off')
            xlabel(vs,                               'Interpreter', 'none', 'FontSize', 18)
            ylabel(param,                            'Interpreter', 'none', 'FontSize', 18)
            title(sprintf('%s vs %s', param, vs),    'Interpreter', 'none', 'FontSize', 20)
            hold off
            
        end
        
        
        % Plot the Stats
        function PlotStat(obj, type, stat, param, vs)
            
%             figure('name', param, 'numbertitle', 'off')
            
            n = numel(type);
            
            p(n) = gobjects;
            
            for ii = 1: n
                
                y_data = obj.Stats.(type{ii}).(param).(stat);
                ind = isnan(y_data);
                y_data(ind) = [];
                
                if nargin > 4
                    x_data = obj.Stats.(type{ii}).(vs).(stat);
                    x_data(ind) = [];
                    
                    ind = min(numel(x_data), numel(y_data));
                    
                    p(ii) = plot(x_data(1:ind), y_data(1:ind));
                    hold on
                else
                    p(ii) = plot(y_data);
                    hold on
                end
                
            end
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting

            
%             legend(p, type, 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 16, 'Box','off')
%             xlabel(vs,                               'Interpreter', 'none', 'FontSize', 18)
%             ylabel(param,                            'Interpreter', 'none', 'FontSize', 18)
%             title(sprintf('%s vs %s', param, vs),    'Interpreter', 'none', 'FontSize', 20)
            hold off
            
        end
        
        
        % Plot the Stats
        function ErrorPlotStat(obj, type, stat, param, vs, errType)
            
            figure('name', param, 'numbertitle', 'off')
            
            n = numel(type);
            
            p(n) = gobjects;
            
            a = 5;
            
            for ii = 1: n
                
                y_data = obj.Stats.(type{ii}).(param).(stat);
                err    = obj.Stats.(type{ii}).(param).(errType);

                ind = isnan(y_data);
                y_data(ind) = [];
               
                if nargin > 4
                    x_data = obj.Stats.(type{ii}).(vs).(stat);
                    x_data(ind) = [];
                    
                    ind = min(numel(x_data), numel(y_data));
                    
                    x_data = x_data(1:15:ind);
                    y_data = y_data(1:15:ind);
                    
                    p(ii) = plot(x_data, y_data);
                    hold on
                    
                    b = numel(y_data);
                
                    in = floor(linspace(a,b,10));
                    x   = x_data(in);
                    y   = y_data(in);
                    err = err(in);
                    errorbar(x,y,err,'.', 'Color', p(ii).Color)
                    
                else
                    p(ii) = plot(y_data);
                    hold on
                end
                
                a = a + 5;
            end
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting

            legend(p, type, 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 16, 'Box','off')
            xlabel(vs,                               'Interpreter', 'none', 'FontSize', 18)
            ylabel(param,                            'Interpreter', 'none', 'FontSize', 18)
            title(sprintf('%s vs %s', param, vs),    'Interpreter', 'none', 'FontSize', 20)
            hold off
            
        end
        
        
    end
end
