function Subplot(x, y1, y2, xl, yl1, yl2, filename, showMarkers)
    % --- Basic settings ---
    lw = 4;  % line width
    fs = 28; % font size
    if nargin < 8, showMarkers = false; end

    % --- Create invisible figure ---
    fig = figure('visible', 'off');

    % --- Marker indices (every 10th point, adjust if needed) ---
    step = max(1, floor(numel(x) / 20));  % about 20 markers total
    m_idx = 1:step:numel(x);

    % ============================================================
    % --- Subplot 1 ---
    subplot(2, 1, 1);
    plot(x, y1, 'LineWidth', lw, 'Color', 'k'); hold on;
    if showMarkers
        plot(x(m_idx), y1(m_idx), 'ob', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
    end
    xlabel(xl, 'FontSize', fs);
    ylabel(yl1, 'FontSize', fs);
    pos = get(gca, 'Position');
    pos(1) = pos(1) + 0.04;
    set(gca, 'Position', pos);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on; grid on;

    % Auto limits
    xmin = min(x(:)); xmax = max(x(:));
    ymin = min(y1(:)); ymax = max(y1(:));
    xpad = 0.05 * (xmax - xmin + eps); ypad = 0.05 * (ymax - ymin + eps);
    xlim([xmin - xpad, xmax + xpad]);
    ylim([ymin - ypad, ymax + ypad]);

    % ============================================================
    % --- Subplot 2 ---
    subplot(2, 1, 2);
    plot(x, y2, 'LineWidth', lw, 'Color', 'k'); hold on;
    if showMarkers
        plot(x(m_idx), y2(m_idx), 'ob', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
    end
    xlabel(xl, 'FontSize', fs);
    ylabel(yl2, 'FontSize', fs);
    pos = get(gca, 'Position');
    pos(1) = pos(1) + 0.04;
    set(gca, 'Position', pos);
    set(gca, 'FontSize', fs, 'LineWidth', lw, 'TickLength', [0.02, 0.02]);
    box on; grid on;

    % Auto limits
    ymin = min(y2(:)); ymax = max(y2(:));
    ypad = 0.05 * (ymax - ymin + eps);
    xlim([xmin - xpad, xmax + xpad]);
    ylim([ymin - ypad, ymax + ypad]);

    % ============================================================
    % --- Save and close ---
    set(gcf, 'Position', [100, 100, 1200, 1000]);
    print(filename, '-dpng', '-r300');
    fprintf('Saved subplot figure to %s\n', filename);
    close(fig);
end
