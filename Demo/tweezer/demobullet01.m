%  DEMOBULLET01 - Tweezer simulation for bullet-shaped particle with correct vertical shape orientation

%  Material properties: water and polystyrene
mat1 = Material(1.33^2, 1);
mat2 = Material(1.59^2, 1);

%  Polar shape parameters of bullet particle
R = 150;   % radius
L = 450;   % length
k = 16;

%  Polar radius equation for bullet
rad = @( t ) ( ...
    ( (R*L) ./ sqrt(R^2 .* cos(t).^2 + L^2 .* sin(t).^2) ).^(-k) + ...
    ( ( (L/R) ./ ( 1e-3 + 0.5 .* (1 - tanh(20 .* cos(t))) .* cos(t)) ).^(-k) ) ...
    ).^(-1/k);

%  Load T-matrix and drag tensor
load bullet.mat
drag = struct('tt', diag(drag.tt).', 'rr', diag(drag.rr).');
disp(drag)

%  Focus lens and beam
NA = 1.3;
lens = optics.lensfocus(mat1, k0, NA, 'nphi', 21, 'ntheta', 20);
e = 2 * normpdf(lens.rho, 0, 1);
e = e(:) * [1, 0, 0]
foc = eval(lens, e);

%  Evaluation points in xy-plane
x = 900 * linspace(-1, 1, 101);
[xx, yy] = ndgrid(x);
e1 = fields(foc, Point(mat1, 1, [xx(:), yy(:), 0*xx(:)]));
e2 = fields(foc, Point(mat1, 1, [xx(:), 0*yy(:), yy(:)]));

%  Plot field in image plane
figure
h1 = subplot(1, 2, 1);
imagesc(x, x, reshape(dot(e1, e1, 2), size(xx)).');
xlabel('x (nm)');
ylabel('y (nm)');
set(gca, 'YDir', 'normal');
axis equal tight
title('XY plane');

h2 = subplot(1, 2, 2);
imagesc(x, x, reshape(dot(e2, e2, 2), size(xx)).');
xlabel('x (nm)');
ylabel('z (nm)');
set(gca, 'YDir', 'normal');
axis equal tight
title('XZ plane');

%  Scatterer setup
fun = @(pos, k0) fields(foc, Point(mat1, 1, pos));
qinc = multipole.incoming(mat1, k0, fun, 'lmax', tmat.lmax, 'diameter', L*2);
disp(tmat.lmax)
scatterer = tweezer.scatterer(tmat, qinc);
fluid = tweezer .fluidpol(drag);

%  Initial position and rotation
pos = [0, 0, 0];
rot = multipole.rotation(105, 'order', 'y');  % zero rotation => pointing up
dt = 3e-5;
nt = 4000;

% Shape initialization

% Define bullet profile along body z-axis
t_shape = linspace(0, pi, 200);    % 0..pi covers one side (tip to tip)
r_shape = rad(t_shape);

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

shape_local = [x_outline; y_outline; z_outline];

% Rotate initial shape and plot it
shape_rot = rot * shape_local;
hold(h1, 'on');
hold(h2, 'on');
h3 = plot(h1, pos(1) + shape_rot(1,:), pos(2) + shape_rot(2,:), 'm-', 'LineWidth', 1.5);
h4 = plot(h2, pos(1) + shape_rot(1,:), pos(3) + shape_rot(3,:), 'm-', 'LineWidth', 1.5);

dirout = zeros(nt, 3);
theta1 = zeros(1, nt);


% === MAIN SIMULATION LOOP ====================================

for it = 1:nt
    [fopt, nopt] = optforce(scatterer, pos, rot);
    [pos, rot] = browniant(fluid, pos, rot, fopt, nopt, dt);

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
    % disp(theta1(it));
    disp(pos)
    drawnow
end

% === FINAL PLOTS =============================================

figure
plot(dirout);
legend('x', 'y', 'z');
xlabel('Time (a.u.)');
ylabel('Orientation direction');
title('Bullet orientation components over time');



% torque_vs_theta.m
thetas = linspace(0, pi, 181);   % 0..180 deg
torque_z = zeros(size(thetas));
pos0 = [0,0,0];
for ii = 1:numel(thetas)
  rot_test = multipole.rotation( thetas(ii)*360/pi, 'order', 'y' );  % convert to deg
  disp(thetas(ii)*180/pi)
  [~, nopt] = optforce(scatterer, pos0, rot_test);
  % project torque onto the axis that rotates the major axis (approx y or x depending)
  % Here we examine torque component about lab x (change if your axis differs)
  torque_z(ii) = nopt(2);    % change index if y is not the rotation axis you want
end

figure;
plot(thetas*360/pi, torque_z);
xlabel('theta (deg)');
ylabel('Torque component');
title('Optical torque vs polar angle (Bullet at lmax = 12)');
grid on;


% Plot cos(theta) vs time
figure
plot(theta1);
xlabel('Time (steps)');
ylabel('cos(\theta)');
titlestr = sprintf('Bullet orientation evolution (L = %.1f, R = %.0f)', L, R);
title(titlestr);
grid on;


