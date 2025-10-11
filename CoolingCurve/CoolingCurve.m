#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window

T_initial = 100;
T_ambient = 25;

k         = 0.05;
dt        = 0.05;
t_end     = 200;

t    = 0 : dt : t_end;
T    = zeros(size(t));
dT   = zeros(size(t));
T(1) = T_initial;

for i = 1:length(t)-1
    %dTdt = -k * (T(i) - T_ambient);
    dTdt = -k * sin(T(i));
    dT(i+1) = dTdt * dt;
    T(i+1) = T(i) + dTdt * dt;
end

% Remove first element to match length
t(1) = [];
T(1) = [];
dT(1) = [];

[ax, h1, h2] = plotyy(t, T, t, dT); 

% Set styles
set(h1, 'LineStyle', '-',  'Color', 'b', 'LineWidth', 2);  % dT (blue)
set(h2, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);  % T (red)

% Label axes
xlabel('Time (s)');
ylabel(ax(1), 'Temperature (°C)');
ylabel(ax(2), 'Temperature Increment (°C)');

% Format axes
set(ax, 'FontSize', 14, 'LineWidth', 1.5, 'TickLength', [0.02 0.02]);

% Save
print("cooling_curve.png", "-dpng");

