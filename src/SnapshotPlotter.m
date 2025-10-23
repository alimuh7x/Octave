classdef SnapshotPlotter < handle
    properties
        filename = '';
        data1 = [];
        data2 = [];
        dual = false;
        snapshot_id = 0;
    end

    methods
        % --- Constructor (filename only) ---
        function obj = SnapshotPlotter(filename)
            obj.filename = filename;
        end

        % --- Add and plot (explicit x, rho, phi, labels) ---
        function add_and_plot(obj, x, field1, varargin)
            % Usage:
            % Single: add_and_plot(x, rho, xlabel, ylabel)
            % Dual:   add_and_plot(x, rho, phi, xlabel, ylabel1, ylabel2)

            obj.snapshot_id = obj.snapshot_id + 1;

            if numel(varargin) == 2
                % Single plot
                xlabel_str  = varargin{1};
                ylabel1_str = varargin{2};
            elseif numel(varargin) == 4
                % Dual subplot
                field2      = varargin{1};
                xlabel_str  = varargin{2};
                ylabel1_str = varargin{3};
                ylabel2_str = varargin{4};
                obj.dual = true;
            else
                error('Invalid number of arguments to add_and_plot');
            end

            % --- Store data ---
            obj.data1 = [obj.data1; field1];
            if obj.dual
                obj.data2 = [obj.data2; field2];
            end

            % --- Plot ---
            fig = figure('visible', 'off'); clf;
            cmap = jet(size(obj.data1,1));

            if obj.dual
                % ======== Subplot 1 (ρ) ========
                subplot(2,1,1); hold on;
                for k = 1:size(obj.data1,1)
                    plot(x, obj.data1(k,:), 'Color', cmap(k,:), 'LineWidth', 1.3);
                end
                xlabel(xlabel_str);
                ylabel(ylabel1_str);
                grid on; box on;
                set(gca, 'FontSize', 14, 'LineWidth', 2, ...
                    'TickLength', [0.02, 0.02]);

                % ======== Subplot 2 (φ) ========
                subplot(2,1,2); hold on;
                for k = 1:size(obj.data2,1)
                    plot(x, obj.data2(k,:), 'Color', cmap(k,:), 'LineWidth', 1.3);
                end
                xlabel(xlabel_str);
                ylabel(ylabel2_str);
                grid on; box on;
                set(gca, 'FontSize', 14, 'LineWidth', 2, ...
                    'TickLength', [0.02, 0.02]);


            else
                % ======== Single plot (ρ) ========
                hold on;
                for k = 1:size(obj.data1,1)
                    plot(x, obj.data1(k,:), 'Color', cmap(k,:), 'LineWidth', 1.3);
                end
                xlabel(xlabel_str);
                ylabel(ylabel1_str);
                grid on; box on;
                set(gca, 'FontSize', 14, 'LineWidth', 2, ...
                    'TickLength', [0.02, 0.02]);
            end

            % --- Save only PNG ---
            print(fig, obj.filename, '-dpng', '-r150');
            close(fig);
        end
    end
end

