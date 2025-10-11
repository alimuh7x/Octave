#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure windows

% ========== Gradient Plot ==========
[x, y] = meshgrid(-5:0.5:5, -5:0.5:5);
T = exp(-0.1*(x.^2 + y.^2));  % corrected exponent

[Tx, Ty] = gradient(T, 0.5, 0.5);  % compute gradients

figure();  % new figure for gradient plot
contourf(x, y, T, 20);  
colorbar;
hold on;
quiver(x, y, Tx, Ty, 'k');  % black arrows
title('2D Gradient Field');
xlabel('x'); ylabel('y');
axis equal tight;
colormap(jet);
print('gradient_plot.png', '-dpng', '-r200');  % save

% ========== Laplacian Plot ==========
[x, y] = meshgrid(-10:0.1:10);       % fine grid for smoothness
r = sqrt(x.^2 + y.^2);
z = sin(r) ./ (r + eps);            % smooth peak-like surface

figure();
surf(x, y, z, 'EdgeColor', 'none'); % remove grid lines
shading interp;                     % smooth color across faces
title('Smooth 3D Surface Plot');
xlabel('x'); ylabel('y'); zlabel('z');

print('Surf.png', '-dpng', '-r200');  % save
