function [ phi, sDL_re, sDL_im ] = RAARupdate_exwv_ReImDsCs_meas( phi, sDL_re, sDL_im, meas, sol )
             
% RAARupdate_exwvmeas_DsCs

%==================================================================================================
% EXITWAVE FROM SYNTHESIS DICTIONARY AND SPARSE CODE

[ Pi_DC_phi_re ] = trainingsetpatches2image( sDL_re.D * sDL_re.C, sDL_re, sol.sz.sz );
[ Pi_DC_phi_im ] = trainingsetpatches2image( sDL_im.D * sDL_im.C, sDL_im, sol.sz.sz );

R_DC_phi = 2 * ( Pi_DC_phi_re + 1i * Pi_DC_phi_im ) - phi;
% R_DC_phi = 2 * ( Pi_DC_phi_re .* exp( 1i * Pi_DC_phi_im )) - phi;

% R_DC_phi = 2 * phi - Pi_DC_phi;



% [ Ip_phi, ~ ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );
% Ip_phi = 2 * sDL.D * sDL.C - Ip_phi;
% 
% [ R_DC_phi ] = trainingsetpatches2image( Ip_phi, sDL, sol.sz.sz );

%==================================================================================================

Pi_m_R_DC_phi = enforce_2DTPAmeas( R_DC_phi, meas, sol.measLPF, sol );

R_m_R_DC_phi = 2 * Pi_m_R_DC_phi - R_DC_phi;

%==================================================================================================

bbeta = 0.5;

Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% Pi_m_phi = enforce_2DTPAmeas( Pi_DC_phi, meas, sol.measLPF, sol );
phi =  bbeta * 0.5 * ( R_m_R_DC_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;