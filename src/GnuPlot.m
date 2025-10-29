classdef GnuPlot < handle
    properties
        filename = 'plot.png';
        plots = {};
        xrange = '';   % e.g. '[0:10]'
        yrange = '';   % e.g. '[0:1]'
        extra_cmds = {};   % user-defined gnuplot commands
    end

    methods
        % --- Constructor ---
        function obj = CurvePlot(filename)
            if nargin >= 1, obj.filename = filename; end
        end

        % --- Add a new curve ---
        function add_plot(obj, x, y, xlabel_text, ylabel_text, label)
            if nargin < 6, label = sprintf('Curve %d', numel(obj.plots) + 1); end
            entry.x = x(:);
            entry.y = y(:);
            entry.xlabel = xlabel_text;
            entry.ylabel = ylabel_text;
            entry.label = label;
            obj.plots{end + 1} = entry;
        end

        % --- Add custom command ---
        function add_command(obj, cmd)
            obj.extra_cmds{end + 1} = cmd;
        end

        % --- Set x-range ---
        function set_xrange(obj, xmin, xmax)
            obj.xrange = sprintf("[%.8g:%.8g]", xmin, xmax);
        end

        % --- Set y-range ---
        function set_yrange(obj, ymin, ymax)
            obj.yrange = sprintf("[%.8g:%.8g]", ymin, ymax);
        end

        % --- Save all plots ---
        function save(obj)
            if isempty(obj.plots)
                error("No plots added. Use add_plot() first.");
            end

            gp = popen("gnuplot", "w");
            first = obj.plots{1};

            % --- Setup ---
            setup = sprintf([ ...
                "set terminal pngcairo size 1200,800 enhanced dashed font 'Verdana,30'\n", ...
                "set output '%s'\n", ...
                "set mxtics 2\n", ...
                "set mytics 2\n", ...
                "set border lw 3\n", ...
                "set style function linespoints\n", ...
                "set style line 1 lw 4 lc rgb 'black'    ps 3 pt 7   pi  50\n", ...
                "set style line 2 lw 4 lc rgb '#DC143C'  ps 3 pt 9   pi  50\n", ...
                "set style line 3 lw 4 lc rgb '#1E90FF'  ps 3 pt 11  pi  50\n", ...
                "set style line 4 lw 4 lc rgb '#008000'  ps 3 pt 13  pi  50\n", ...
                "set style line 5 lw 4 lc rgb '#0072BD'  ps 3 pt 5   pi  50\n", ...
                "set style line 6 lw 4 lc rgb '#D95319'  ps 3 pt 15  pi  50\n", ...
                "set style line 7 lw 4 lc rgb '#77AC30'  ps 3 pt 19  pi  50\n", ...
                "set style line 8 lw 4 lc rgb '#DC143C'  ps 3 pt 23  pi  50 \n", ...
                "set style line 9 lw 4 lc rgb '#1E90FF'  ps 3 pt 29  pi  50 \n", ...
                "set style line 10 lw 4 lc rgb '#008000' ps 3 pt 31  pi  50 \n", ...
                "set style line 11 lw 4 lc rgb '#FFD700' ps 3 pt 62  pi  50 \n", ...
                "set style line 12 lw 4 lc rgb '#FF69B4' ps 3 pt 28  pi  50 \n", ...
                "set style line 13 lw 4 lc rgb '#8A2BE2' ps 3 pt 7   pi  50 \n", ...
                "set style line 14 lw 4 lc rgb '#00CED1' ps 3 pt 9   pi  50 \n", ...
                "set style line 14 lw 4 lc rgb '#FF4500' ps 3 pt 13  pi  50 \n", ...
                "set style line 15 lw 4 lc rgb '#9400D3' ps 3 pt 31  pi  50 \n", ...
                "set xlabel '%s'\n", ...
                "set ylabel '%s'\n", ...
                "set key bottom\n", ...
                "set grid\n" ...
            ], obj.filename, first.xlabel, first.ylabel);

            % --- Add optional ranges ---
            if ~isempty(obj.xrange)
                setup = [setup, sprintf("set xrange %s\n", obj.xrange)];
            end
            if ~isempty(obj.yrange)
                setup = [setup, sprintf("set yrange %s\n", obj.yrange)];
            end

            % --- Add user commands ---
            for i = 1:numel(obj.extra_cmds)
                setup = [setup, obj.extra_cmds{i}, "\n"];
            end

            fprintf(gp, "%s", setup);
            fflush(gp);

            % --- Plot command ---
            fprintf(gp, "plot ");
            for k = 1:numel(obj.plots)
                plt = obj.plots{k};
                ls = mod(k-1,15)+1;
                if k < numel(obj.plots)
                    fprintf(gp, "'-' using 1:2 with lines ls %d title '%s',\\\n", ls, plt.label);
                else
                    fprintf(gp, "'-' using 1:2 with lines ls %d title '%s'\n", ls, plt.label);
                end
            end
            fflush(gp);

            % --- Stream data + flush ---
            for k = 1:numel(obj.plots)
                plt = obj.plots{k};
                data = [plt.x(:), plt.y(:)];
                fprintf(gp, "%.8f %.8f\n", data');
                fprintf(gp, "e\n");
                fflush(gp);
            end

            fprintf(gp, "unset output\n");
            fflush(gp);
            pclose(gp);

            fprintf("✅ Saved %d curve(s) to %s\n", numel(obj.plots), obj.filename);
        end
    end
end
