function Subplot(x, y1, y2, xl, yl1, yl2, filename)
    % --- Basic settings ---
    lw = 4;              % line width
    fs = 24;             % font size

    % --- Create invisible figure ---
    fig = figure('visible', 'off');

    % ============================================================
    % --- Subplot 1 ---
    subplot(2, 1, 1);
    plot(x, y1, 'LineWidth', lw);
    xlabel(xl, 'FontSize', fs);
    ylabel(yl1, 'FontSize', fs);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on; grid on;

    % Auto limits with padding
    xmin = min(x(:)); xmax = max(x(:));
    ymin = min(y1(:)); ymax = max(y1(:));
    xrange = xmax - xmin; yrange = ymax - ymin;
    if xrange == 0, xrange = 1; end
    if yrange == 0, yrange = 1; end
    xpad = 0.05 * xrange; ypad = 0.05 * yrange;
    xlim([xmin - xpad, xmax + xpad]);
    ylim([ymin - ypad, ymax + ypad]);

    % ============================================================
    % --- Subplot 2 ---
    subplot(2, 1, 2);
    plot(x, y2, 'LineWidth', lw);
    xlabel(xl, 'FontSize', fs);
    ylabel(yl2, 'FontSize', fs);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on; grid on;

    % Auto limits with padding
    ymin = min(y2(:)); ymax = max(y2(:));
    yrange = ymax - ymin;
    if yrange == 0, yrange = 1; end
    ypad = 0.05 * yrange;
    xlim([xmin - xpad, xmax + xpad]);
    ylim([ymin - ypad, ymax + ypad]);

    % ============================================================
    % --- Save and close ---
    set(gcf, 'Position', [100, 100, 1200, 1000]);  % wide layout
    print(filename, '-dpng', '-r300');
    fprintf('Saved subplot figure to %s\n', filename);
    close(fig);
end

