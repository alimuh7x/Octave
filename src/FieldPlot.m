classdef FieldPlot < handle
    % ---------------------------------------------------------------------
    % FieldPlot: minimal plotting helper for scalar + vector fields
    % Automatically normalizes arrows for cleaner visualization.
    %
    % Example:
    %   p = FieldPlot("potential.png");
    %   p.add_Image(x, y, phi, "Electric potential φ", ColorMaps.turbo);
    %   p.add_Arrows(X, Y, Ex, Ey, 15, 0.07, Colors.black);
    %   p.save();
    % ---------------------------------------------------------------------

    properties
        filename
        fig_handle
    end

    methods
        % --- Constructor ---
        function obj = FieldPlot(filename)
            if nargin < 1
                error("Usage: p = FieldPlot('filename.png')");
            end
            obj.filename = filename;
            obj.fig_handle = figure('visible', 'off');
            hold on;
        end

        % --- Add scalar image ---
        function add_Image(obj, x, y, field, Title, cmap)
            if nargin < 6, cmap = ColorMaps.jet; end
            if nargin < 5, Title = ""; end

            imagesc(x, y, field);
            set(gca, "YDir", "normal");
            axis image; colorbar; colormap(cmap);

            if ~isempty(Title)
                title(Title);
            end
        end

        % --- Add electric field arrows (auto-normalized) ---
        % Usage: add_Arrows(X, Y, U, V, skip, scale, color)
        function add_Arrows(obj, X, Y, U, V, skip, scale, color)
            if nargin < 8, color = Colors.black; end
            if nargin < 7, scale = 0.05; end
            if nargin < 6, skip = 6; end

            % Sample points for clarity
            Xs = X(1:skip:end, 1:skip:end);
            Ys = Y(1:skip:end, 1:skip:end);
            Us = U(1:skip:end, 1:skip:end);
            Vs = V(1:skip:end, 1:skip:end);

            % --- Normalize vectors ---
            mag = sqrt(Us.^2 + Vs.^2);
            mag(mag == 0) = 1;  % avoid division by zero
            Us = Us ./ mag;
            Vs = Vs ./ mag;

            % Plot arrows (normalized, scaled)
            quiver(Xs, Ys, scale * Us, scale * Vs, 0, color, "LineWidth", 5);
        end

        % --- Save figure ---
        function save(obj)
            print(obj.filename, "-dpng");
            hold off;
            close(obj.fig_handle);
            fprintf("Saved figure to %s\n", obj.filename);
        end
    end
end
