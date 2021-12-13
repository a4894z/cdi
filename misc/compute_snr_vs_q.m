function [ snr_reconst, phot_px ] = compute_snr_vs_q( meas_modulus )

sz = size( meas_modulus );

meas_region = ( meas_modulus ~= 0 );

%====================================================================================================================================================

%keep this for the time being, compare speeds for defining it here vs passing it
[ c, r ] = meshgrid( 1 : sz( 2 ), 1 : sz( 1 ) );

%========

%the q = 0 value:
jj = 1; kk = 1;

o_diam_circ = ( ( ( sz( 2 ) / 2 + 1 - c ) ./ jj ) .^ 2 + ( ( sz( 1 ) / 2 + 1 - r ) ./ jj ) .^2 < 1 );

snr_reconst( kk ) = sum( sum( o_diam_circ .* meas_modulus)) / ( sum( sum( sqrt( o_diam_circ .* meas_modulus ))));
phot_px( kk )     = sum( sum( o_diam_circ .* meas_modulus));

%========

jj_max = round(sz( 2 )/2); 
annuls_smplng = 1;

%after q = 0, which corresponds to jj = 1, we now need to loop from:
loop_range = ((1 + annuls_smplng) : annuls_smplng : jj_max);

%========

snr_reconst(2 : (1 + length(loop_range))) = 0;
phot_px(2 : (1 + length(loop_range))) = 0;

%{
snr_num(2 : (1 + length(loop_range))) = 0;
snr_denom(2 : (1 + length(loop_range))) = 0;
%}

%========

kk = 2;

tmpc = ( sz( 2 ) / 2 + 1 - c );
tmpr = ( sz( 1 ) / 2 + 1 - r );

for jj = loop_range
  
  o_diam_circ = ( ( tmpc ./ jj ) .^ 2        + ( tmpr ./ jj )         .^ 2 < 1 );
  i_diam_circ = ( ( tmpc ./ ( jj - 1 )) .^ 2 + ( tmpr ./ ( jj - 1 ) ) .^ 2 < 1 );
  
  temp1 = meas_region .* ( o_diam_circ - i_diam_circ );
  temp2 = temp1 .* meas_modulus;
  temp3 = sum( sum( temp2 ));
   
  phot_px( kk )     = temp3 / sum( sum( temp1 ));
  snr_reconst( kk ) = temp3 / ( 1e-7 + sum( sum( sqrt( temp2 ))));

  kk = kk + 1;
  
end

