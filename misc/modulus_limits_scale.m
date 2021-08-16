function [ V_bounded ] = modulus_limits_scale( V, lim )

%================================================================
% rescale so that array magnitude is between 0 and 1
%================================================================

abs_V = abs( V );

min_abs = min( abs_V( : ));
abs_V = abs_V - min_abs; 

max_abs = max( abs_V( : ));
abs_V = abs_V / ( 1e-7 + max_abs );

%================================================================
% rescale so that we're between our desired magnitude range
%================================================================

V_bounded = ( abs_V * ( lim( 2 ) - lim( 1 )) + lim( 1 )) .* exp( 1i * angle( V ));
