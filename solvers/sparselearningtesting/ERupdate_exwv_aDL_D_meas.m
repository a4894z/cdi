function [ phi, aDL ] = ERupdate_exwv_aDL_D_meas( phi, aDL, meas, sol )

%==================================================================================================
% EXIT WAVE MEASUREMENT CONSTRAINT UPDATE:

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

% tmp0 = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% phi = 2 * tmp0 - phi;

% phi = phi .* sol.probe.S;

%==================================================================================================
% ANALYSIS SPARSE DICTIONARY UPDATE:
 
[ Ip, aDL ] = image2trainingsetpatches( phi, aDL, sol.sz );

%=======================

% aDL.D = aDL.C * pinv( Ip );    % minimizer for min_Da || Da * Ip - Ca ||^2_F

[ u, s, v ] = svd( aDL.C * pinv( Ip ));
s = eye( sol.aDL.Na, sol.aDL.Ni );
aDL.D = u * s * v';

% % rescale so that WHOLE MATRIX has max value of one
% aDL.D = aDL.D / max( abs( aDL.D(:) ) );   

% % rescale so that each analysis dictionary row has max value of one
% aDL.D = aDL.D ./ ( 1e-7 + repmat( max( abs( aDL.D ), [], 2), [ 1, sol.aDL.Ni ] ));

%=======================

[ phi ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sol.sz );

%==================================================================================================
