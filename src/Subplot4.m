function Subplot4(array, names, filename, showMarkers)
    % --- Basic settings ---
    lw = 1;  % line width
    fs = 10; % font size
    if nargin < 4, showMarkers = false; end

    % --- Extract x and y ---
    x = array(:, 1);
    yArray = array(:, 2:end);
    yNames = names(2:end);
    xName  = names(1);

    % --- Create invisible figure ---
    fig = figure('visible', 'off', 'Renderer', 'painters');

    nPlots = size(yArray, 2);
    for i = 1:nPlots
        subplot(2, 2, i);

        % --- Remove invalid data ---
        mask = isfinite(x) & isfinite(yArray(:, i));
        xf = x(mask);
        yf = yArray(mask, i);

        % --- Marker indices ---
        step = max(1, floor(numel(xf) / 20));
        m_idx = 1:step:numel(xf);

        % --- Plot line + optional markers ---
        plot(xf, yf, 'k', 'LineWidth', lw); hold on;
        if showMarkers
            plot(xf(m_idx), yf(m_idx), 'ob', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
        end

        % --- Adjust subplot spacing ---
        set(gca, 'Position', get(gca, 'Position') + [0.0, 0.0, 0.02, 0.02]);
        if (i == 2 || i == 4)
            set(gca, 'Position', get(gca, 'Position') + [0.03, 0.0, 0.0, 0.0]);
        end

        % --- Labels ---
        xlabel(xName, 'FontSize', fs, 'Interpreter', 'tex');
        ylabel(yNames(i), 'FontSize', fs, 'Interpreter', 'tex');

        % --- Styling ---
        set(gca, 'FontSize', fs, 'LineWidth', 1, 'TickLength', [0.02, 0.02]);
        box on;

        % --- Limits with padding ---
        xmin = min(xf(:)); xmax = max(xf(:));
        ymin = min(yf(:)); ymax = max(yf(:));
        xpad = 0.05 * (xmax - xmin + eps);
        ypad = 0.05 * (ymax - ymin + eps);
        xlim([xmin - xpad, xmax + xpad]);
        ylim([ymin - ypad, ymax + ypad]);
    end

    % --- Save ---
    print(filename, '-dpng', '-r300');
    fprintf('Saved 2x2 subplot figure to %s\n', filename);
    close(fig);
end
