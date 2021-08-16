function [ V_bounded ] = modulus_limits_project( V, lim )

abs_x = abs( V );

temp0 = ( abs_x >= lim( 1 ) );
temp1 = ( abs_x <= lim( 2 ) );
temp2 = temp0 .* temp1;

V_bounded = temp2 .* V + ( lim( 2 ) * not( temp1 ) + lim( 1 ) * not( temp0 ) ) .* exp( 1i * angle( V ));
