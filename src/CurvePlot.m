classdef CurvePlot < handle
    properties
        filename = '';
        fig;
        sub_count = 0;
        total_sub = 0;
        legends = {};
        line_width = 2;
        font_size = 15;
        xmin = inf;
        xmax = -inf;
        ymin = inf;
        ymax = -inf;
    end

    methods
        % --- Constructor ---
        function obj = CurvePlot(filename)
            if nargin > 0
                obj.filename = filename;
            end

            % Force thick lines globally
            set(0, 'defaultlinelinewidth', obj.line_width);

            % Create hidden figure
            obj.fig = figure('visible', 'off');
        end

        % --- Declare total subplots before plotting ---
        function set_total_sub(obj, n)
            obj.total_sub = n;
        end

        function update_limits(obj, x, y)
            obj.xmin = min(obj.xmin, min(x(:)));
            obj.xmax = max(obj.xmax, max(x(:)));
            obj.ymin = min(obj.ymin, min(y(:)));
            obj.ymax = max(obj.ymax, max(y(:)));
        end


        % --- Plot a single curve ---
        function plot_curve(obj, x, y, xl, yl, varargin)
            if isempty(findall(obj.fig, 'Type', 'axes'))
                subplot(1,1,1);
            end
            plot(x, y, varargin{:});
            obj.update_limits(x, y);
            xlabel(xl); ylabel(yl);
            set(gca, 'FontSize', obj.font_size, 'TickLength', [0.02, 0.02]);
            set(gca, 'linewidth', obj.line_width);
            set(gca, 'clipping', 'on');
            hold on;

            obj.legends{end+1} = yl;
        end

        % --- Add another curve to same plot ---
        function add_plot(obj, x, y, xl, yl, varargin)
            plot(x, y, varargin{:});
            obj.update_limits(x, y);
            xlabel(xl); ylabel(yl);
            set(gca, 'FontSize', obj.font_size, 'TickLength', [0.02, 0.02]);
            set(gca, 'linewidth', obj.line_width);
            set(gca, 'clipping', 'on');
            hold on;
            obj.legends{end+1} = yl;
        end

        % --- Add subplots correctly ---
        function add_sub(obj, x, y, xl, yl, varargin)
            if obj.total_sub == 0
                error('Use set_total_sub(n) before adding subplots.');
            end
            obj.sub_count = obj.sub_count + 1;
            subplot(obj.total_sub, 1, obj.sub_count);
            plot(x, y, varargin{:});
            obj.update_limits(x, y);
            xlabel(xl); ylabel(yl);
            set(gca, 'FontSize', obj.font_size, 'TickLength', [0.02, 0.02]);
            set(gca, 'linewidth', obj.line_width);
            set(gca, 'clipping', 'on');
        end

        % --- Save figure ---
        function save(obj)
            if obj.sub_count == 0 && ~isempty(obj.legends)
                legend(obj.legends, 'fontsize', 16, 'box', 'off', 'location', 'northeast');
            end
            if isfinite(obj.xmin)
                % Compute data ranges
                xrange = obj.xmax - obj.xmin;
                yrange = obj.ymax - obj.ymin;
            
                % Add 5% margin on each side (adjust if needed)
                xpad = 0.05 * xrange;
                ypad = 0.05 * yrange;
            
                % Handle constant data (avoid zero range)
                if xrange == 0, xpad = 1; end
                if yrange == 0, ypad = 1; end
            
                % Apply limits
                xlim([obj.xmin - xpad, obj.xmax + xpad]);
                ylim([obj.ymin - ypad, obj.ymax + ypad]);
            end
            box on;  % redraw axes border



            if ~isempty(obj.filename)
                set(gca, 'Position', [0.15 0.15 0.75 0.75]);  % manual tightening if needed
                print(obj.filename, '-dpng', '-r300', '-tight');  % use -tight
                fprintf('Saved to %s\n', obj.filename);
                close(obj.fig);
            end
        end
    end
end

