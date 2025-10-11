#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window

% Grid setup
N = 200;
x = linspace(0, 10, N);
dx = x(2) - x(1);
dt = 0.001;
steps = 1000;

% Parameters
D = 0.1;                 % diffusion coefficient
f = sin(pi * x / 10).^2;  % initial field, non-negative
f = f / trapz(x, f);      % normalize mass

% Parabolic potential
V = 0.5 * (x - 5).^2;
dVdx = gradient(V, dx);   % numerical gradient

% Storage for plotting
figure;
colors = jet(6);
plot_interval = floor(steps / 5);
plot(x, f, 'color', colors(1,:), 'linewidth', 2); hold on;

plot_id = 2;

for step = 1:steps
    % Compute gradients
    df_dx  = gradient(f, dx);
    d2f_dx2 = gradient(df_dx, dx);

    diffusion = D * d2f_dx2; % diffusion: D ∂^2f/∂x^2 
    adv = gradient(f .* dVdx, dx);     % advection: ∂(f ∂V/∂x)/∂x

    Totalrate = diffusion + adv; % net rate of change
    
    % Update
    f = f + dt * (Totalrate);
    f = max(f, 0);       % keep non-negative
    f = f / trapz(x, f); % conserve total mass

    % Plot at intervals
    if mod(step, plot_interval) == 0
        plot(x, f, 'color', colors(plot_id,:), 'linewidth', 1.5);
        plot_id += 1;
    end
end

% Plot potential
plot(x, V / max(V), 'k--', 'linewidth', 2);
legend('Initial', '', '', '', '', 'Normalized V(x)');
xlabel('x'); ylabel('f(x)');
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);

% Save figure
print("diffusion_advection.png", "-dpng");

