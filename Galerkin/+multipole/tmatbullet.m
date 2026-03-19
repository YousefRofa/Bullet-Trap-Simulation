function t = tmatbullet( mat1, mat2, R, L, k0, varargin )

%  set up solver for particles with symmetry of revolution
rad = @(t) R + ((1 + cos(t))/2) .* ...  
    ((L.*R)./sqrt(R.^2.*cos(t).^2 + L.^2.*sin(t).^2) - R);
pol = multipole.polsolver( mat1, mat2, rad, varargin{ : } );
%  evaluate T-matrices
t = arrayfun( @( k0 ) eval( pol, k0 ), k0, 'uniform', 1 );
