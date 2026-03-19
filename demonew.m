%  DEMOFORCE02 - Force on nanobullet.

%  material properties
mat1 = Material( 1.33 ^ 2, 1 );
mat2 = Material( 1.59 ^ 2, 1 );
%  material vector
mat = [ mat1, mat2 ];

%  polar shape of particle
scale = 4;
R = scale * 150;
L = scale * 450;
rad = @( t ) R + ( ( 1 + cos( t ) ) / 2 ) .* ...  
    ( ( L * R ) ./ sqrt( R ^ 2 * cos( t ) .^ 2 + L ^ 2  * sin( t ) .^ 2 ) - R );

%  discretized boundary for nanosphere
u = linspace( 0, 2 * pi, 50 );
t = pi * linspace( 0, 1, 70 ) .^ 2;
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
%  BEM solver 
bem = fill( galerkin.bemsolver( tau, 'order', [] ), k0 );

%  focus lens
NA = 1.0;
lens = optics.lensfocus( mat1, k0, NA, 'nphi', 21, 'ntheta', 20 );
%  incoming fields
e = normpdf( lens.rho, 0, 1 );
e = e( : ) * [ 1, 0, 0 ];
%  planewave decomposition of focal fields
foc = eval( lens, e );
foc2 = optics.rotx( 180 ) * foc;

%  particle positions
zout = 4e3 * linspace( -1, 1, 201 );
pos = zout( : ) * [ 0, 0, 1 ];
%  optical force and torque from BEM solution
fun = @( i ) foc( tau, 'shift', pos( i, : ) );
q = arrayfun( fun, 1 : size( pos, 1 ), 'uniform', 1 );
q = struct( 'tau', tau, 'k0', k0, 'e', horzcat( q.e ), 'h', horzcat( q.h ) );
sol1 = solve( bem, q );
[ fopt1, nopt1 ] = optforce( sol1 );

%  optical force and torque from BEM solution
fun = @( i ) foc2( tau, 'shift', pos( i, : ) );
q = arrayfun( fun, 1 : size( pos, 1 ), 'uniform', 1 );
q = struct( 'tau', tau, 'k0', k0, 'e', horzcat( q.e ), 'h', horzcat( q.h ) );
sol2 = solve( bem, q );
[ fopt2, nopt2 ] = optforce( sol2 );
%  flip forces to account for different propagation direction of fields
fopt2 = - flip( fopt2, 1 );


%  final plot

% figure

plot( 1e-3 * zout, fopt1( :, 3 ), 1e-3 * zout, fopt2( :, 3 ) );  hold on
set( gca, 'ColorOrderIndex', 1 );

xlabel( 'z (\mum)' );
ylabel( 'Optical force (pN)' );
