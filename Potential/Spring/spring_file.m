%#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "on");  % hide figure window
clc; clear; clearvars;

N = 100;                     % number of points
x = linspace(0, 10, N);
dx = x(2) - x(1);
f = sin(pi * x / 10);        % initial field

k = 1;                       % spring constant
M = 1;                       % mobility
V = 0.5 * (x - 5).^4;        % external potential
dt = 0.001;                  % time step
steps = 1000;
plot_interval = 200;         % plot every 200 steps

% Plot initial field
figure;
colors = jet(steps/plot_interval + 1);

plot(x, f, 'color', colors(1,:), 'linewidth', 2); hold on;

plot_count = 1;

for step = 1:steps
    f_new = f;
    for i = 2:N-1
        lap = (f(i+1) - 2*f(i) + f(i-1)) / dx^2;
        f_new(i) = f(i) + dt * M * (k * lap - 2*V(i)*f(i));
    end
    f = f_new;

    if mod(step, plot_interval) == 0
        plot_count += 1;
        plot(x, f, 'color', colors(plot_count,:), 'linewidth', 1.5);
        drawnow;
    end
end

% Plot normalized potential
plot(x, V / max(V), 'k--', 'linewidth', 2);


xlabel('x'); ylabel('f(x)');
legend('Initial f(x)', '', '', 'Normalized V(x)');
axis([0 10 -0.2 1.2]);
set   (gca, 'FontSize', 30);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);

% Save final multi-curve plot
print("spring_chain_snapshots.png", "-dpng");
