function [ V ] = make_3Dgaussian( N, mu, sdev )

[ c, r, p ] = meshgrid( 1 : N( 2 ), 1 : N( 1 ), 1 : N( 3 ) );

p = single( p ); 
c = single( c ); 
r = single( r ); 

V_p = exp( -1 * ( p - mu( 3 ) ).^2 / ( 1e-7 + 2 * sdev( 3 ) )^2 );
V_c = exp( -1 * ( c - mu( 2 ) ).^2 / ( 1e-7 + 2 * sdev( 2 ) )^2 );
V_r = exp( -1 * ( r - mu( 1 ) ).^2 / ( 1e-7 + 2 * sdev( 1 ) )^2 );

V = V_p .* V_r .* V_c;

