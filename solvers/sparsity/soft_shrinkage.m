function [ B, W ] = soft_shrinkage( A, W, lam )

%==================================================================================================

if ( lam < 0 )

  warning( 'must use positive real thresholding parameter, doing nothing' ); 
  return
  
end

%==================================================================================================

% determine where ( abs_A - lam ) > 0
W = ( W > lam );
% W = ( abs_A > lam );

% remember phase
% phs_A = ( A ./ ( 1e-7 + abs_A ) );
phs_A = angle( A );

%================================================

% soft thresholding
abs_A = abs( A );
B = ( abs_A - lam ) .* W .* exp( 1i * phs_A );  

% % soft thresholding (might be faster):
% B = ( A - lam * phs_A ) .* W; 

%==================================================================================================