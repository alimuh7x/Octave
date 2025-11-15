function SinglePlot(x, y, xl, yl, filename, showMarkers)
    % --- Basic settings ---
    lw = 2;              % line width
    fs = 15;             % font size
    if nargin < 6, showMarkers = false; end

    step = max(1, floor(numel(x) / 20));  % about 20 markers total
    m_idx = 1:step:numel(x);
    % --- Create figure ---
    fig = figure('visible', 'off');
    plot(x, y, 'LineWidth', lw, 'Color', 'k'); hold on;
    if showMarkers
        plot(x(m_idx), y(m_idx), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    end
    xlabel(xl, 'FontSize', fs);
    ylabel(yl, 'FontSize', fs);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on;

    % --- Auto-adjust limits with small padding ---
    xmin = min(x(:)); xmax = max(x(:));
    ymin = min(y(:)); ymax = max(y(:));
    xrange = xmax - xmin; yrange = ymax - ymin;
    if xrange == 0, xrange = 1; end
    if yrange == 0, yrange = 1; end
    xpad = 0.05 * xrange; ypad = 0.05 * yrange;
    xlim([xmin - xpad, xmax + xpad]);
    ylim([ymin - ypad, ymax + ypad]);

    % --- Save plot ---
    set(gca, 'Position', [0.20 0.15 0.75 0.75]);  % manual tightening if needed
    print(filename, '-dpng');
    fprintf('Saved plot to %s\n', filename);
    close(fig);
end

