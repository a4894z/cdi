function [ V ] = make_3Dspheroid( N, diam )

[ c, r, p ] = meshgrid( -0.5 * N(2) + 1 : 0.5 * N(2) - 0, ...
                        -0.5 * N(1) + 1 : 0.5 * N(1) - 0, ...
                        -0.5 * N(3) + 1 : 0.5 * N(3) - 0 );
                    
r = single( r );                  
c = single( c ); 
p = single( p ); 

V = (( r ./ ( 0.5 * diam(1) )).^2 + ( c ./ ( 0.5 * diam(2) )).^2 + ( p ./ ( 0.5 * diam(3) )).^2 ) < 1;

V = single( V );

end
