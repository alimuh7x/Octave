#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window

disp("✅ Code Runner is working");
x = linspace(0, 10, 1000);
V = 0.5 * (x - 5).^4;
f = sin(2* pi * x / 10);

dt = 0.005;
nSteps = 2000;
nPlots = 5;
plot_every = floor(nSteps / nPlots);

figure;
colors = jet(nPlots);  % color map for multiple lines
plot_count = 0;

m = -0.001 * V .* f;
plot(x, m, 'r', 'linewidth', 1.5);  hold on;
plot(x, V / max(V), 'k--', 'linewidth', 2);
xlabel('x'); ylabel('- V * f');

set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
% Save final multi-curve plot
print("V_f.png", "-dpng");
% clear previous plots
clf;


plot(x, f, 'k', 'linewidth', 1.5); hold on;

for step = 1:nSteps
    f = f - dt * V .* f;
    %f = f / norm(f);

    if mod(step, plot_every) == 0
        plot_count += 1;
        plot(x, f, 'color', colors(plot_count, :), 'linewidth', 1.5); hold on;
    end
end

disp("✅ Plotting f(x) under parabolic potential");
% Plot potential for reference
% legend_entries = arrayfun(@(k) sprintf("Step %d", k * plot_every), 1:nPlots, 'UniformOutput', false);
%legend([legend_entries, {"Normalized V(x)"}]);
title('Evolution of f(x) under parabolic potential');
xlabel('x'); ylabel('f(x)');
axis([0 10 -1.2 1.2]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);

% Save final multi-curve plot
print("f_evolution_snapshots.png", "-dpng");

disp("✅ Finished");
