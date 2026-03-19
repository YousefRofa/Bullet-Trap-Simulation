%  DEMOBULLET_BEM - Tweezer simulation for bullet-shaped particle using BEM (no T-matrix)

clear; clc; close all;

% =====================================================
% 1. MATERIAL AND SHAPE DEFINITION
% =====================================================

%  Material properties: water and polystyrene
mat1 = Material(1.33^2, 1);
mat2 = Material(1.59^2, 1);
mat = [mat1, mat2];

%  Polar shape of particle
scale = 4;
R = scale * 150;      % nm
L = scale * 450;      % nm
rad = @(t) R + ((1 + cos(t)) / 2) .* ...
    ((L * R) ./ sqrt(R^2 * cos(t).^2 + L^2 * sin(t).^2) - R);

% =====================================================
% 2. GEOMETRY AND BEM SETUP
% =====================================================

u = linspace(0, 2*pi, 50);
t = pi * linspace(0, 1, 70).^2;
p = trispheresegment(u, t, 2);

[~, t] = cart2sph(p.verts(:,1), p.verts(:,2), p.verts(:,3));
t = 0.5*pi - t;
p.verts = rad(t) .* p.verts;

% Boundary elements with linear shape functions
tau = BoundaryEdge(mat, p, [2,1]);

% Wavenumber of light in vacuum
lambda0 = 520;              % nm
k0 = 2 * pi / lambda0;

% BEM solver
bem = fill(galerkin.bemsolver(tau, 'order', []), k0);

% =====================================================
% 3. FOCUS AND BEAM
% =====================================================

NA = 1.0;
lens = optics.lensfocus(mat1, k0, NA, 'nphi', 21, 'ntheta', 20);
e = normpdf(lens.rho, 0, 1);
e = e(:) * [1, 0, 0];
foc = eval(lens, e);
foc2 = optics.rotx(180) * foc;

% =====================================================
% 4. PARTICLE POSITIONS AND BEM FORCE CALCULATION
% =====================================================

zout = 4e3 * linspace(-1, 1, 201);
pos = zout(:) * [0, 0, 1];

fun = @(i) foc(tau, 'shift', pos(i,:));
q = arrayfun(fun, 1:size(pos,1), 'uniform', 1);
q = struct('tau', tau, 'k0', k0, 'e', horzcat(q.e), 'h', horzcat(q.h));
sol1 = solve(bem, q);
[fopt1, nopt1] = optforce(sol1);

fun = @(i) foc2(tau, 'shift', pos(i,:));
q = arrayfun(fun, 1:size(pos,1), 'uniform', 1);
q = struct('tau', tau, 'k0', k0, 'e', horzcat(q.e), 'h', horzcat(q.h));
sol2 = solve(bem, q);
[fopt2, nopt2] = optforce(sol2);
fopt2 = -flip(fopt2, 1);

figure;
plot(1e-3*zout, fopt1(:,3), 1e-3*zout, fopt2(:,3), 'LineWidth', 1.5);
xlabel('z (\mum)');
ylabel('Optical force (pN)');
legend('Forward', 'Backward');
title('Optical Force along z');
grid on;

% =====================================================
% 5. SIMPLE DYNAMICS SIMULATION
% =====================================================
t_shape = linspace(0, pi, 200);    % 0..pi covers one side (tip to tip)
r_shape = rad(t_shape);

%  Initial position and rotation
pos = [0, 0, 0];
rot = multipole.rotation(270, 'order', 'y');  % zero rotation => pointing up
dt = 3e-5;
nt = 4000;

% In local coordinates:
%   z is along bullet axis
%   r is radial distance from z-axis
x_shape = r_shape .* sin(t_shape);
y_shape = zeros(size(x_shape));
z_shape = r_shape .* cos(t_shape);   % so shape points along +Z

% Mirror to make full revolution (2D outline around z-axis)
x_outline = [x_shape, fliplr(-x_shape)];
y_outline = [y_shape, fliplr(y_shape)];
z_outline = [z_shape, fliplr(z_shape)];

shape_local = 0.2*[x_outline; y_outline; z_outline];

% Rotate initial shape and plot it
shape_rot = rot * shape_local;
hold(h1, 'on');
hold(h2, 'on');
h3 = plot(h1, pos(1) + shape_rot(1,:), pos(2) + shape_rot(2,:), 'm-', 'LineWidth', 1.5);
h4 = plot(h2, pos(1) + shape_rot(1,:), pos(3) + shape_rot(3,:), 'm-', 'LineWidth', 1.5);

dirout = zeros(nt, 3);
theta1 = zeros(1, nt);

for it = 1:nt
    % Simplified optical force (approx near center)
    fopt = [0; 0; mean(fopt1(:,3)) * 1e-12];  % Convert pN to N
    nopt = [0; 0; 0];  % neglect torque for simplicity

    [pos, rot] = browniant_custom(pos, rot, fopt, nopt, dt, drag_tt, drag_rr);
 % Rotate bullet shape according to current orientation
    shape_rot = rot * shape_local;

    % Update shape position in both subplots
    set(h3, 'XData', pos(1) + shape_rot(1,:), 'YData', pos(2) + shape_rot(2,:));
    set(h4, 'XData', pos(1) + shape_rot(1,:), 'YData', pos(3) + shape_rot(3,:));

    % Store orientation
    dirout(it, :) = rot * [0; 0; 1];

    % Compute and log cos(theta)
    z_body = [0; 0; 1];
    z_lab  = rot * z_body;
    theta1(it) = z_lab(3);
    disp(theta1(it));

    drawnow
end


figure
plot(dirout);
legend('x', 'y', 'z');
xlabel('Time (a.u.)');
ylabel('Orientation direction');
title('Bullet orientation components over time');

% Plot cos(theta) vs time
figure
plot(theta1);
xlabel('Time (steps)');
ylabel('cos(\theta)');
titlestr = sprintf('Bullet orientation evolution (L = %.1f, R = %.0f)', L, R);
title(titlestr);
grid on;

% =====================================================
% 6. CUSTOM BROWNIAN INTEGRATOR
% =====================================================

function [pos, rot] = browniant_custom(pos, rot, fopt, nopt, dt, drag_tt, drag_rr)
    % Brownian step with translational and rotational diffusion
    kB = 1.380649e-23;   % Boltzmann constant (J/K)
    T = 300;             % Temperature (K)

    M_t = inv(drag_tt);
    M_r = inv(drag_rr);

    % Deterministic drift
    v = M_t * fopt;
    omega = M_r * nopt;

    % Stochastic terms
    xi_t = sqrtm(2 * kB * T * M_t * dt) * randn(3,1);
    xi_r = sqrtm(2 * kB * T * M_r * dt) * randn(3,1);

    pos = pos + (v + xi_t) * dt;
    rot = rot * axang2rotm([xi_r'/norm(xi_r + 1e-12), norm(xi_r)]);  % small rotation
end
