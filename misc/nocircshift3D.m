function [ Xshift ] = nocircshift3D( X, shift_vec_rc, fillreg )

if ~exist( 'fillreg', 'var' ), fillreg = 0; end

% non-circular shift an array:

Xs = padarray( X, abs( shift_vec_rc ), fillreg );
Xs = circshift( Xs, shift_vec_rc );

sz  = size( X ); 
sz2 = size( Xs );

csr = sz2(1)/2 + 1 - sz(1)/2 : sz2(1)/2 + sz(1)/2;
csc = sz2(2)/2 + 1 - sz(2)/2 : sz2(2)/2 + sz(2)/2;
csd = sz2(3)/2 + 1 - sz(3)/2 : sz2(3)/2 + sz(3)/2;

Xshift = Xs( csr, csc, csd );
    
