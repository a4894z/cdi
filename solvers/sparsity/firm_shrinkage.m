function [ B, W1 ] = firm_shrinkage( A, W, lam1, lam2 )

% firm shrinkage is:
            
%                { 0                                                      if      | x | <= lam1
% Pi_{fs}[ x ] = { sgn( x ) lam2 ( | x | - lam1 ) / ( lam2 - lam1 )       if      lam1 < | x | <= lam2
%                { x                                                      if      | x | > lam2
%
% It reduces to hard shrinkage for lam2 = lam1, and soft shrinkage for lam2 --> inf


%==================================================================================================

if ( lam1 > lam2 )

  warning( 'must use lam1 > lam2 for firm shrinkage, doing nothing' ); 
  return
  
end


if ( lam1 < 0 ) || ( lam2 < 0 ) 

  warning( 'must use positive real thresholding parameters for firm shrinkage, doing nothing' ); 
  return
  
end

%==================================================================================================

abs_A = abs( A );

W2 = ( W > lam2 );
W1 = ( W <= lam1 );
W12 = not( W2 ) & not( W1 );

%B = W12 .* ( ( lam2 / ( lam2 - lam1 ) ) * ( abs_A - lam1 ) ) .* ( A ./ ( 1e-7 + abs_A ) ) + A .* W2;
B = W12 .* ( ( lam2 / ( lam2 - lam1 ) ) * ( abs_A - lam1 ) ) .* exp( 1i * angle( A )) + A .* W2;

% DO I HAVE TO RECOMPUTE THIS, OR DO I ALREADY HAVE IT ABOVE?
W1 = ( B ~= 0 );

%==================================================================================================
