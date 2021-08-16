function [ phi, sDL ] = RAARupdate_exwv_DsCs_meas( phi, sDL, meas, sol )
             
% RAARupdate_exwvmeas_DsCs

%==================================================================================================
% EXITWAVE FROM SYNTHESIS DICTIONARY AND SPARSE CODE

[ Pi_DC_phi ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sol.sz.sz );

R_DC_phi = 2 * Pi_DC_phi - phi;
% R_DC_phi = 2 * phi - Pi_DC_phi;



% [ Ip_phi, ~ ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );
% Ip_phi = 2 * sDL.D * sDL.C - Ip_phi;
% 
% [ R_DC_phi ] = trainingsetpatches2image( Ip_phi, sDL, sol.sz.sz );

%==================================================================================================

Pi_m_R_DC_phi = enforce_2DTPAmeas( R_DC_phi, meas, sol.measLPF, sol );

R_m_R_DC_phi = 2 * Pi_m_R_DC_phi - R_DC_phi;

%==================================================================================================

bbeta = 0.7;

Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% Pi_m_phi = enforce_2DTPAmeas( Pi_DC_phi, meas, sol.measLPF, sol );
phi =  bbeta * 0.5 * ( R_m_R_DC_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;