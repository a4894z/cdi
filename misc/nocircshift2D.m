function [ Xshift ] = nocircshift2D( X, shift_vec_rc )

% non-circular shift an array:

% Xs = zeropadarray( X, abs( shift_vec_rc ));
% Xs = padarray( X( :, :, 3 ), abs( shift_vec_rc ));
Xs = padarray( X, abs( shift_vec_rc ));

Xs = circshift( Xs, shift_vec_rc );

sz  = size( X ); 
sz2 = size( Xs );

csr = ( sz2(1)/2 + 1 - sz(1)/2 ) : ( sz2(1)/2 + sz(1)/2 );
csc = ( sz2(2)/2 + 1 - sz(2)/2 ) : ( sz2(2)/2 + sz(2)/2 );

Xshift = Xs( csr, csc, : );
    
