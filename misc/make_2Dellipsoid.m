function [ V ] = make_2Dellipsoid( N, diam )

[ c, r ] = meshgrid( -0.5 * N( 2 ) : 0.5 * N( 2 ) - 1,  -0.5 * N( 1 ) : 0.5 * N( 1 ) - 1 );
c = single( c ); 
r = single( r ); 

%create ellipsoid:
V = ( ( c ./ ( 0.5 * diam( 2 ) )).^2 + ( r ./ ( 0.5 * diam( 1 ) )).^2 ) < 1;
V = single( V );

end
