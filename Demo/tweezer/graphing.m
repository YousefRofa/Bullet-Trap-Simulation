% Polar Plot of Custom Radius Function
% ------------------------------------

clear; clc; close all;

% Parameters
R = 200;       % Example value (change as needed)
L = 600;       % Example value (change as needed)
k = 16;      % Given parameter

% Define the radius function (as a function handle)
rad = @(t) ( ...
    ( (R .* L) ./ sqrt(R^2 .* cos(t).^2 + L^2 .* sin(t).^2) ).^(-k) + ...
    ( ( (L./R) ./ ( 1e-3 + 0.5 .* (1 - tanh(20 .* cos(t))) .* abs(cos(t))) ).^(-k) ) ...
    ).^(-1/k);

% Define angle range (0 to 2*pi)
t = linspace(0, 2*pi, 1000);

% Compute radius values
r = rad(t);

% Create polar plot
figure;
polarplot(t, r, 'LineWidth', 2);
title('Polar Plot of Custom Radius Function');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
grid on;

% Optional: also show Cartesian equivalent
figure;
plot(r .* cos(t), r .* sin(t), 'LineWidth', 2);
axis equal;
xlabel('x');
ylabel('y');
title('Cartesian View of Polar Curve');
grid on;
