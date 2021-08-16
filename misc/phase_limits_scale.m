function [ V_bounded ] = phase_limits_scale( V, lim )

%================================================================
% rescale so that phase is between 0 and 1
%================================================================

angle_V = angle( V );

min_phs = min( angle_V( : ));
angle_V = angle_V - min_phs; 

max_phs = max( angle_V( : ));
angle_V = angle_V / ( 1e-7 + max_phs );

%================================================================
% rescale so that we're between our desired phase range
%================================================================

angle_V = angle_V * ( lim( 2 ) - lim( 1 )) + lim( 1 );

V_bounded = abs( V ) .* exp( 1i * angle_V );







