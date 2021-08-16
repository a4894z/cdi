function [ f2, f1, W ] = trainingsetpatches2image( Y, DL, sz )

Y1 = reshape( Y, [ DL.wy, DL.wx, DL.Np ] );
% Y1 = Y1 - repmat( mean( mean( Y1 )), [ DL.wy, DL.wx ] );

if DL.mean0 == true && isfield( DL , 'theta' )
    
    Y1 = Y1 - repmat( mean( mean( Y1 )), [ DL.wy, DL.wx ] );
    Y1 = Y1 + reshape( repmat( DL.theta, [ DL.Ni, 1 ] ), [ DL.wy, DL.wx, DL.Np ] );
    
end

W = zeros( sz( 1 ), sz( 2 ) );
f1 = zeros( sz( 1 ), sz( 2 ) );

for ii = 1 : DL.Np
    
    x = DL.Xp( :, :, ii ); 
    y = DL.Yp( :, :, ii );
    
    slc = y + ( x - 1 ) * sz( 1 );
    f1( slc ) = f1( slc ) + Y1( :, :, ii );
    W( slc ) = W( slc ) + 1;
    
end

f2 = f1 ./ W;

f2( isnan( f2 )) = 0;
f2( isinf( f2 )) = 0;