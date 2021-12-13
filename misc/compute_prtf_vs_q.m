function [ prtf ] = compute_prtf_vs_q( meas_modulus, soln_modulus )

%====================================================================================================================================================

meas_region = ( meas_modulus ~= 0 );

sz = size( soln_modulus );

[ xx, yy ] = meshgrid( 1 : sz( 2 ), 1 : sz( 1 ) );

%====================================================================================================================================================
% q = 0 value:
%=============

jj = 1;
o_diam_circ = ( (( 0.5 * sz( 2 ) + 1 - xx ) ./ jj ) .^ 2 + ( ( 0.5 * sz( 1 ) + 1 - yy ) ./ jj ) .^ 2 < 1 );

kk = 1;
prtf( kk ) = sum( sum( o_diam_circ .* soln_modulus )) / ( 1e-7 + sum( sum( o_diam_circ .* meas_modulus )));

kk = kk + 1;

%====================================================================
% after q = 0, which corresponds to jj = 1, we now need to loop from:
%====================================================================

skip = 1;

loop_range = ( ( 1 + skip ) : skip : round( sz( 2 ) / 2 ) );

prtf( 2 : ( 1 + length( loop_range ))) = 0;

tmpx = ( sz( 2 ) / 2 + 1 - xx );
tmpy = ( sz( 1 ) / 2 + 1 - yy );
  
for jj = loop_range
  
  o_diam_circ = ( ( tmpx ./ jj )        .^2 + ( tmpy ./ jj )        .^2 < 1 );
  i_diam_circ = ( ( tmpx ./ ( jj - 1 )) .^2 + ( tmpy ./ ( jj - 1 )) .^2 < 1 );
  
  temp1 = meas_region .* ( o_diam_circ - i_diam_circ );
  
  prtf(kk) = sum( sum( temp1 .* soln_modulus )) / ( 1e-7 + sum( sum( temp1 .* meas_modulus )) );

  kk = kk + 1;

end

%====================================================================================================================================================
