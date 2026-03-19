%  DEMOTMAT01 - T-matrix for dielectric nanobullet and single wavelength.

%  material properties
mat1 = Material( 1.33 ^ 2, 1 );
mat2 = Material( 1.59 ^ 2, 1 );
%  material vector
mat = [ mat1, mat2 ];

%  polar shape of particle
R = 150;
L = 450;
k = 16;

% PARTICLE SHAPE
% Ensure that the shape is the same as that in demobullet1

% Perfect Elipse 
% rad = @( t )  ( L * R ) ./ sqrt( R ^ 2 * cos( t ) .^ 2 + L ^ 2  * sin( t ) .^ 2 );


% Approximated Shape
% rad = @( t ) R + ( ( 1 + cos( t ) ) / 2 ) .* ...  
%     ( ( L * R ) ./ sqrt( R ^ 2 * cos( t ) .^ 2 + L ^ 2  * sin( t ) .^ 2 ) - R );

% More accurate approximation of shape, including hyperbolic fucntions
% Note that this shape approximation is not to be fully trusted sinc
% hyperbolic funciton may miss up the maths for laser reflection

rad = @( t ) ( ...
    ( (R*L) ./ sqrt(R^2 .* cos(t).^2 + L^2 .* sin(t).^2) ).^(-k) + ...
    ( ( (L/R) ./ ( 1e-3 + 0.5 .* (1 - tanh(20 .* cos(t))) .* cos(t)) ).^(-k) ) ...
    ).^(-1/k);

%  discretized boundary for nanosphere
u = linspace( 0, 2 * pi, 30 );
t = pi * linspace( 0, 1, 30 ) .^ 2;
p = trispheresegment( u, t, 2 );
%  scale radius
[ ~, t ] = cart2sph( p.verts( :, 1 ), p.verts( :, 2 ), p.verts( :, 3 ) );
t = 0.5 * pi - t;
p.verts = rad( t ) .* p.verts;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 2, 1 ] );     
%  plot( tau, 'EdgeColor', 'b' );

%  wavenumber of light in vacuum
k0 = 2 * pi / 520;
%  BEM solver and T-matrix solver
lmax = 10;
bem = fill( galerkin.bemsolver( tau, 'order', [] ), k0 );
tsolver = multipole.tsolver( mat, 1, lmax );
%  T-matrix
sol = bem \ tsolver( tau, k0 );
tmat = eval( tsolver, sol );

% %  additional information for H5 file
% info = multipole.h5info( tau );
% info.name = "Bullet";
% info.description = "Single bullet and single wavelength";
% info.matname = [ "Water", "Polystyrene" ];
% %  save T-matrix
% fout = 'tmatrix_bullet.h5';
% h5save( tmat, fout, info );

drag = tweezer.stokesdrag( p, 'quad', triquad( 3 ) );

%  save T-matrix and drag tensor
save bullet k0 tmat drag 
