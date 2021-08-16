function [ X ] = pinv_sparse( A )

[ U, S, V ] = svds( A, max(size( A )));

s = diag( S );
 
if nargin < 2 
    tol = max( size( A )) * eps( norm( s, inf ));
end

r1 = sum( s > tol ) + 1;

V( :, r1 : end ) = [];
U( :, r1 : end ) = [];
s( r1 : end ) = [];
s = 1 ./ s( : );
X = ( V .* s.' ) * U';