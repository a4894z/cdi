function [ Y, DL ] = image2trainingsetpatches( f, DL, sz )

Y = [];

%==================================================================================================

% % w x w patch indices:
% [ dX, dY ] = meshgrid( 0 : ( DL.w - 1), 0 : ( DL.w - 1 ));

% wx x wy patch indices:
[ dX, dY ] = meshgrid( 0 : ( DL.wx - 1 ), 0 : ( DL.wy - 1 ));

%==================================================================================================

if strcmp( DL.patch, 'random' )

    % random patch locations:
    
    % use whole image:
    x = floor( rand( 1, 1, DL.Np ) * ( sz( 2 ) - DL.wx ) ) + 1;
    y = floor( rand( 1, 1, DL.Np ) * ( sz( 1 ) - DL.wy ) ) + 1;

%     % subregion of image:
%     x = floor( 64 + rand( 1, 1, DL.Np ) * ( 128 - DL.wx ) ) + 1;
%     y = floor( 64 + rand( 1, 1, DL.Np ) * ( 128 - DL.wy ) ) + 1;

    % Extract lots of patches y_j \in R^n, and store them in a matrix Y = (y_j)_{j=1}^m.
%     DL.Xp = repmat( dX, [ 1, 1, DL.Np ] ) + repmat( x, [ DL.w, DL.w, 1 ] );
%     DL.Yp = repmat( dY, [ 1, 1, DL.Np ] ) + repmat( y, [ DL.w, DL.w, 1 ] );
    DL.Xp = repmat( dX, [ 1, 1, DL.Np ] ) + repmat( x, [ DL.wy, DL.wx, 1 ] );
    DL.Yp = repmat( dY, [ 1, 1, DL.Np ] ) + repmat( y, [ DL.wy, DL.wx, 1 ] );
    
elseif strcmp( DL.patch, 'raster' )

    % define regularly space positions for the extraction of patches:

    % q > 0 is an overlap factor ( so that setting q = w implies no overlap ). 
    qx = round( DL.wx * ( 1 - DL.ol ) );
    qy = round( DL.wy * ( 1 - DL.ol ) );

    [ x, y ] = meshgrid( round( 1 : qx : ( ( sz( 2 ) - DL.wx / 2 ) - 0 ) ), round( 1 : qy : ( ( sz( 1 ) - DL.wy / 2 ) - 0 )));
%     [ x, y ] = meshgrid( round( 1 : DL.q : ( ( sz( 2 ) - DL.w / 2 ) - 0 ) ), round( 1 : DL.q : ( ( sz( 1 ) - DL.w / 2 ) - 0 )));

    DL.Np = size( x( : ), 1 );

%     DL.Xp = repmat( dX, [ 1, 1, DL.Np ] ) + repmat( reshape( x( : ), [ 1, 1, DL.Np ] ), [ DL.w, DL.w, 1 ] );
%     DL.Yp = repmat( dY, [ 1, 1, DL.Np ] ) + repmat( reshape( y( : ), [ 1, 1, DL.Np ] ), [ DL.w, DL.w, 1 ] );

    DL.Xp = repmat( dX, [ 1, 1, DL.Np ] ) + repmat( reshape( x( : ), [ 1, 1, DL.Np ] ), [ DL.wy, DL.wx, 1 ] );
    DL.Yp = repmat( dY, [ 1, 1, DL.Np ] ) + repmat( reshape( y( : ), [ 1, 1, DL.Np ] ), [ DL.wy, DL.wx, 1 ] );

    % ensure boundary conditions (reflection):  
    % ASH: what about trying wrap around, like in circshift?
    tmp0 = ( DL.Yp > sz( 1 ) );
    if any( tmp0( : ))
        DL.Yp( DL.Yp > sz( 1 ) ) = 2 * sz( 1 ) - DL.Yp( DL.Yp > sz( 1 ) );
    end

    tmp0 = ( DL.Xp > sz( 2 ) );
    if any( tmp0( : ))
        DL.Xp( DL.Xp > sz( 2 ) ) = 2 * sz( 2 ) - DL.Xp( DL.Xp > sz( 2 ) );
    end

end

%==================================================================================================

if ~isempty( f )
    
    % Extract the Np patches Y. This notation implicitly column stacks the image, uses each
    % column of the argument to select pixels, and returns the original dim of image.
    Y = f( DL.Yp + ( DL.Xp - 1 ) * sz( 1 ) );

    % 2d image array --> column stack vector
    Y = reshape( Y, [ DL.Ni, DL.Np ] );

    %============================================
    % only keep those with largest energy.

    % [ tmp, I ] = sort( sum( Y .^ 2 ), 'descend' );
    % Y = Y( :, I( 1 : DL.Np ));

    %============================================
    % remove the mean for each patch

    if DL.mean0 == true
        
        DL.theta = mean( Y );
        Y = Y - repmat( DL.theta, [ DL.Ni, 1 ] );

    end

end

%==================================================================================================
