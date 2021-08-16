function [ the_rectangle ] = make_3Dcuboid( N, len )

% create rectangle:
slcx = round( 0.5 * N( 2 ) + 0 - 0.5 * len( 2 ) : 0.5 * N( 2 ) + 0.5 * len( 2 ) );
slcy = round( 0.5 * N( 1 ) + 0 - 0.5 * len( 1 ) : 0.5 * N( 1 ) + 0.5 * len( 1 ) );
slcz = round( 0.5 * N( 3 ) + 0 - 0.5 * len( 3 ) : 0.5 * N( 3 ) + 0.5 * len( 3 ) );

% if the side lengths are too big for the array size, just set the size of the 3d cuboid to be the same as the array size
if any( slcx < 0 ), slcx = 1 : N( 2 ); end
if any( slcy < 0 ), slcy = 1 : N( 1 ); end
if any( slcz < 0 ), slcz = 1 : N( 3 ); end

the_rectangle =  zeros( N( 1 ), N( 2 ), N( 3 ), 'single' );
the_rectangle( slcy, slcx, slcz ) = 1;

end
