function [ the_rectangle ] = make_rectangle( N, len )

% create rectangle:
slcx = round( 0.5 * N( 2 ) + 1 - 0.5 * len( 2 ) : 0.5 * N(2) + 0.5 * len( 2 ) );
slcy = round( 0.5 * N( 1 ) + 1 - 0.5 * len( 1 ) : 0.5 * N(1) + 0.5 * len( 1 ) );

% the side lengths are too big for the array size, so just set the size of the rectangle to be the
% same as the array size
if any( slcx < 0 ), slcx = 1 : N(2); end
if any( slcy < 0 ), slcy = 1 : N(1); end

the_rectangle = single( zeros( N( 1 ), N( 2 ) ));
the_rectangle( slcy, slcx ) = 1;

end
