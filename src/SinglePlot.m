function SinglePlot(x, y, xl, yl, filename)
    % --- Basic settings ---
    lw = 2;              % line width
    fs = 15;             % font size

    % --- Create figure ---
    fig = figure('visible', 'off');
    plot(x, y, 'LineWidth', lw);
    xlabel(xl, 'FontSize', fs);
    ylabel(yl, 'FontSize', fs);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on; grid on;

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
    print(filename, '-dpng', '-r300');
    fprintf('Saved plot to %s\n', filename);
    close(fig);
end

