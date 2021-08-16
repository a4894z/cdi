function [ phi, sDL ] = RAARupdate_exwv_sDL_C_meas( phi, sDL, meas, sol )

% RAARupdate_exwvmeas_sDLsparsecode

%==================================================================================================
% SYNTHESIS SPARSE CODE UPDATE

[ Ip, sDL ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );

sDL.C = sDL_sparsecode_update( Ip, sDL );
% sDL.D = sDL_dictionary_update( Ip, sDL );

[ Pi_C_phi ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sol.sz.sz );

R_C_phi = 2 * Pi_C_phi - phi;

%==================================================================================================

Pi_m_R_C_phi = enforce_2DTPAmeas( R_C_phi, meas, sol.measLPF, sol );

R_m_R_C_phi = 2 * Pi_m_R_C_phi - R_C_phi;

%==================================================================================================

bbeta = 0.5;

Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% Pi_m_phi = enforce_2DTPAmeas( R_C_phi, meas, sol.measLPF, sol );
phi =  bbeta * 0.5 * ( R_m_R_C_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;



