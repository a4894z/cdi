function [ Vs, phi ] = RAARupdate_exwv_sparseFDxy_meas_V2( Vs, sFDxy, meas, sol )


%==================================================================================================

abs_Vs.y = abs( Vs.x );
abs_Vs.x = abs( Vs.y );

sFDxy.sparselvl = round( sol.sz.rc * 0.12 );

tmp0 = sqrt( abs_Vs.y .^ 2 + abs_Vs.x .^ 2 );
thresh = find_thresh_from_sparsitylevel( tmp0, sFDxy.sparselvl );

[ Pi_s_phi.y, Wy ] = hard_shrinkage( Vs.y, tmp0, thresh );
[ Pi_s_phi.x, Wx ] = hard_shrinkage( Vs.x, tmp0, thresh );
% [ Pi_s_phi.y, Wy ] = soft_shrinkage( Vs.y, tmp0, thresh );
% [ Pi_s_phi.x, Wx ] = soft_shrinkage( Vs.x, tmp0, thresh );




% 
% abs_Vs.y = abs( Vs.x );
% abs_Vs.x = abs( Vs.y );
% 
% sFDxy.sparselvl = round( sol.sz.rc * 0.09 );
%     
% threshr = find_thresh_from_sparsitylevel( abs_Vs.y, sFDxy.sparselvl );
% threshc = find_thresh_from_sparsitylevel( abs_Vs.x, sFDxy.sparselvl );
% 
% [ Pi_s_phi.y, Wy ] = hard_shrinkage( Vs.y, abs_Vs.y, threshr );
% [ Pi_s_phi.x, Wx ] = hard_shrinkage( Vs.x, abs_Vs.x, threshc );
% % [ Pi_s_phi.y, Wy ] = soft_shrinkage( Vs.y, abs_Vs.y, threshr );
% % [ Pi_s_phi.x, Wx ] = soft_shrinkage( Vs.x, abs_Vs.x, threshc );
%  






%     Pi_s_re_phi.y = real( Vs.y );
%     Pi_s_re_phi.x = real( Vs.x );
%     
%     Pi_s_im_phi.y = imag( Vs.y );
%     Pi_s_im_phi.x = imag( Vs.x );
% 
%     sFDxy.sparselvl = round( sol.sz.rc * 0.09 );
%     
%     tmpre = sqrt( Pi_s_re_phi.y .^ 2 + Pi_s_re_phi.x .^ 2 );
%     tmpim = sqrt( Pi_s_im_phi.y .^ 2 + Pi_s_im_phi.x .^ 2 );
% 
%     threshre = find_thresh_from_sparsitylevel( tmpre, sFDxy.sparselvl );
%     threshim = find_thresh_from_sparsitylevel( tmpim, sFDxy.sparselvl );
% 
%     [ Pi_s_re_phi.y, Wrey ] = hard_shrinkage( Pi_s_re_phi.y, tmpre, threshre );
%     [ Pi_s_re_phi.x, Wrex ] = hard_shrinkage( Pi_s_re_phi.x, tmpre, threshre );
%     [ Pi_s_im_phi.y, Wimy ] = hard_shrinkage( Pi_s_im_phi.y, tmpim, threshim );
%     [ Pi_s_im_phi.x, Wimx ] = hard_shrinkage( Pi_s_im_phi.x, tmpim, threshim );
%     
%     Pi_s_phi.y = Pi_s_re_phi.y + 1i * Pi_s_im_phi.y;
%     Pi_s_phi.x = Pi_s_re_phi.x + 1i * Pi_s_im_phi.x;

   



        

%==================

R_s_Vs.x = 2 * Pi_s_phi.x - Vs.x;
R_s_Vs.y = 2 * Pi_s_phi.y - Vs.y;

%==================================================================================================
    
R_s_Vs_phi = iedgedetect_FDxy( R_s_Vs, sFDxy, sol.sz );

Pi_m_R_s_Vs_phi = enforce_2DTPAmeas( R_s_Vs_phi, meas, sol.measLPF, sol );

[ Pi_m_R_s_Vs ] = edgedetect_FDxy( Pi_m_R_s_Vs_phi, sFDxy );

%==================

R_m_R_s_Vs.x = 2 * Pi_m_R_s_Vs.x - R_s_Vs.x;
R_m_R_s_Vs.y = 2 * Pi_m_R_s_Vs.y - R_s_Vs.y;

%==================

[ phiVs ] = iedgedetect_FDxy( Vs, sFDxy, sol.sz );

phi = enforce_2DTPAmeas( phiVs, meas, sol.measLPF, sol );

[ Pi_m_Vs ] = edgedetect_FDxy( phi, sFDxy );

%==================================================================================================

bbeta = 0.5;

Vs.x = 0.5 * ( 1 - bbeta ) * ( R_m_R_s_Vs.x + Vs.x ) + bbeta * Pi_m_Vs.x;
Vs.y = 0.5 * ( 1 - bbeta ) * ( R_m_R_s_Vs.y + Vs.y ) + bbeta * Pi_m_Vs.y;

%==================================================================================================

