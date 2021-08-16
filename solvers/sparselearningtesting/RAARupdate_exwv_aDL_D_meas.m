function [ phi, aDL ] = RAARupdate_exwv_aDL_D_meas( phi, aDL, meas, sol )
              

% 0.5 * bbeta * ( R_m[ R_s[ phi ]] + phi ) + ( 1 - bbeta ) * Pi_m[ phi ]


%==================================================================================================
% ORIGINAL IMAGE TO IMAGE PATCHES

[ Ip, aDL ] = image2trainingsetpatches( phi, aDL, sol.sz );

%==================================================================================================
% ANALYSIS SPARSE DICTIONARY UPDATE:

% aDL.D = aDL.C * pinv( Ip );    % minimizer for min_Da || Da * Ip - Ca ||^2_F

[ u, s, v ] = svd( aDL.C * pinv( Ip ));
s = eye( sol.aDL.Na, sol.aDL.Ni );
aDL.D = u * s * v';

%==================================================================================================

% % rescale so that WHOLE MATRIX has max value of one
% aDL.D = aDL.D / max( abs( aDL.D(:) ) );   

%=======================

% % rescale so that each analysis dictionary row has max value of one
% aDL.D = aDL.D ./ ( 1e-7 + repmat( max( abs( aDL.D ), [], 2), [ 1, sol.aDL.Ni ] ));

%==================================================================================================

[ Pi_D_phi ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sol.sz );

R_D_phi = 2 * Pi_D_phi - phi;

%==================================================================================================

Pi_m_R_D_phi = enforce_2DTPAmeas( R_D_phi, meas, sol.measLPF, sol );

R_m_R_D_phi = 2 * Pi_m_R_D_phi - R_D_phi;

%==================================================================================================

Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

bbeta = 0.5;

phi = 0.5 * bbeta * ( R_m_R_D_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;






