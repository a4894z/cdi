function [ phi, sDL ] = RAARupdate_Ds_Cs( phi, sDL, meas, sol )
             
%==================================================================================================
% SYNTHESIS SPARSE DICTIONARY UPDATE:

[ Ip, sDL ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );

[ sDL.D ] = sDL_dictionary_update( Ip, sDL.C, sDL.D );

[ Pi_D_phi ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sol.sz.sz );

R_D_phi = 2 * Pi_D_phi - phi;

%==================================================================================================

Pi_m_R_D_phi = enforce_2DTPAmeas( R_D_phi, meas, sol.measLPF, sol );

R_m_R_D_phi = 2 * Pi_m_R_D_phi - R_D_phi;

%==================================================================================================

bbeta = 0.7;

% Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
Pi_m_phi = enforce_2DTPAmeas( Pi_D_phi, meas, sol.measLPF, sol );
phi = 0.5 * bbeta * ( R_m_R_D_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;


