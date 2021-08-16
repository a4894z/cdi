function [ phi ] = RAARupdate_exwv_sparseFDxy_meas( phi, sFDxy, meas, sol )


Vs = edgedetect_FDxy( phi, sFDxy );





% abs_Vs.y = abs( Vs.x );
% abs_Vs.x = abs( Vs.y );
% 
% sFDxy.sparselvl = round( sol.sz.rc * 0.09 );
% 
% tmp0 = sqrt( abs_Vs.y .^ 2 + abs_Vs.x .^ 2 );
% thresh = find_thresh_from_sparsitylevel( tmp0, sFDxy.sparselvl );
% 
% [ Vs.y, Wy ] = hard_shrinkage( Vs.y, tmp0, thresh );
% [ Vs.x, Wx ] = hard_shrinkage( Vs.x, tmp0, thresh );
% % [ Vs.y, Wy ] = soft_shrinkage( Vs.y, tmp0, thresh );
% % [ Vs.x, Wx ] = soft_shrinkage( Vs.x, tmp0, thresh );






% 
% abs_Vs.y = abs( Vs.x );
% abs_Vs.x = abs( Vs.y );
% 
% sFDxy.sparselvl = round( sol.sz.rc * 0.06 );
%     
% threshr = find_thresh_from_sparsitylevel( abs_Vs.y, sFDxy.sparselvl );
% threshc = find_thresh_from_sparsitylevel( abs_Vs.x, sFDxy.sparselvl );
% 
% [ Vs.y, Wy ] = hard_shrinkage( Vs.y, abs_Vs.y, threshr );
% [ Vs.x, Wx ] = hard_shrinkage( Vs.x, abs_Vs.x, threshc );
% % [ Vs.y, Wy ] = soft_shrinkage( Vs.y, abs_Vs.y, threshr );
% % [ Vs.x, Wx ] = soft_shrinkage( Vs.x, abs_Vs.x, threshc );
    








    re_Vs.y = real( Vs.y );
    re_Vs.x = real( Vs.x );
    
    im_Vs.y = imag( Vs.y );
    im_Vs.x = imag( Vs.x );

    
    sFDxy.sparselvl = round( sol.sz.rc * 0.09 );
    
    
    tmpre = sqrt( re_Vs.y .^ 2 + re_Vs.x .^ 2 );
    tmpim = sqrt( im_Vs.y .^ 2 + im_Vs.x .^ 2 );

    threshre = find_thresh_from_sparsitylevel( tmpre, sFDxy.sparselvl );
    threshim = find_thresh_from_sparsitylevel( tmpim, sFDxy.sparselvl );

    [ re_Vs.y, Wrey ] = hard_shrinkage( re_Vs.y, tmpre, threshre );
    [ re_Vs.x, Wrex ] = hard_shrinkage( re_Vs.x, tmpre, threshre );
    [ im_Vs.y, Wimy ] = hard_shrinkage( im_Vs.y, tmpim, threshim );
    [ im_Vs.x, Wimx ] = hard_shrinkage( im_Vs.x, tmpim, threshim );
    
    Vs.y = re_Vs.y + 1i * im_Vs.y;
    Vs.x = re_Vs.x + 1i * im_Vs.x;














[ Pi_s_phi ] = iedgedetect_FDxy( Vs, sFDxy, sol.sz );
        
        


R_s_phi = 2 * Pi_s_phi - phi;

%==================================================================================================
    
  
Pi_m_R_s_phi = enforce_2DTPAmeas( R_s_phi, meas, sol.measLPF, sol );

Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% Pi_m_phi = enforce_2DTPAmeas( Pi_s_phi, meas, sol.measLPF, sol );
%     

R_m_R_s_phi = 2 * Pi_m_R_s_phi - R_s_phi;

%==================================================================================================

aalpha = 0.666;
bbeta = 0.666;

phi = 0.5 * aalpha * bbeta * ( R_m_R_s_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi + ( 1 - aalpha )  * Pi_s_phi;
% phi = 0.5 * ( 1 - aalpha ) * ( 1 - bbeta ) * ( R_m_R_s_phi + phi ) + bbeta * Pi_m_phi;



phi = lpf_gauss( phi, 0.3 * sol.sz.sz );

%==================================================================================================

