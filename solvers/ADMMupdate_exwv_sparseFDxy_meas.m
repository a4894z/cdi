function [ phi, Vm, Vs, lam ] = ADMMupdate_exwv_sparseFDxy_meas( Vm, lam, sFDxy, meas, sol )

    Vs.x = Vm.x - lam.x;
    Vs.y = Vm.y - lam.y;

    %===================

%     abs_Vs.y = abs( Vs.x );
%     abs_Vs.x = abs( Vs.y );
% 
%     
%     sFDxy.sparselvl = round( sol.sz.rc * 0.3 );
%     
%     tmp0 = sqrt( abs_Vs.y .^ 2 + abs_Vs.x .^ 2 );
%     thresh = find_thresh_from_sparsitylevel( tmp0, sFDxy.sparselvl );
% 
%     [ Vs.y, Wy ] = hard_shrinkage( Vs.y, tmp0, thresh );
%     [ Vs.x, Wx ] = hard_shrinkage( Vs.x, tmp0, thresh );
% %     [ Vs.y, Wy ] = soft_shrinkage( Vs.y, tmp0, thresh );
% %     [ Vs.x, Wx ] = soft_shrinkage( Vs.x, tmp0, thresh );
%     
% 
% %     
% % % %     threshr = find_thresh_from_sparsitylevel( abs_Vs.y, sFDxy.sparselvl );
% % % %     threshc = find_thresh_from_sparsitylevel( abs_Vs.x, sFDxy.sparselvl );
% % % % 
% % % % %     [ Vs.y, Wy ] = hard_shrinkage( Vs.y, abs_Vs.y, threshr );
% % % % %     [ Vs.x, Wx ] = hard_shrinkage( Vs.x, abs_Vs.x, threshc );
% % % %     [ Vs.y, Wy ] = soft_shrinkage( Vs.y, abs_Vs.y, threshr );
% % % %     [ Vs.x, Wx ] = soft_shrinkage( Vs.x, abs_Vs.x, threshc );
    









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
     



    

%     re_Vs.y = real( Vs.y );
%     re_Vs.x = real( Vs.x );
%     
%     im_Vs.y = imag( Vs.y );
%     im_Vs.x = imag( Vs.x );
% 
%     abs_re_Vs.y = abs( re_Vs.y );
%     abs_re_Vs.x = abs( re_Vs.x );
%     abs_im_Vs.y = abs( im_Vs.y );
%     abs_im_Vs.x = abs( im_Vs.x );
%     
%     sFDxy.sparselvl = round( sol.sz.rc * 0.10 );
%       
%     threshrey = find_thresh_from_sparsitylevel( abs_re_Vs.y, sFDxy.sparselvl );
%     threshrex = find_thresh_from_sparsitylevel( abs_re_Vs.x, sFDxy.sparselvl );
%     threshimy = find_thresh_from_sparsitylevel( abs_im_Vs.y, sFDxy.sparselvl );
%     threshimx = find_thresh_from_sparsitylevel( abs_im_Vs.x, sFDxy.sparselvl );
% 
%     [ re_Vs.y, Wrey ] = hard_shrinkage( re_Vs.y, abs_re_Vs.y, threshrey );
%     [ re_Vs.x, Wrex ] = hard_shrinkage( re_Vs.x, abs_re_Vs.x, threshrex );
%     [ im_Vs.y, Wimy ] = hard_shrinkage( im_Vs.y, abs_im_Vs.y, threshimy );
%     [ im_Vs.x, Wimx ] = hard_shrinkage( im_Vs.x, abs_im_Vs.x, threshimx );
%     
%     Vs.y = re_Vs.y + 1i * im_Vs.y;
%     Vs.x = re_Vs.x + 1i * im_Vs.x;
    
    




    
    %===================
    
    Vm.x = Vs.x + lam.x;
    Vm.y = Vs.y + lam.y;
    
    %===================
    
    [ phi ] = iedgedetect_FDxy( Vm, sFDxy, sol.sz );

    phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
    
%     phi = lpf_gauss( phi, 0.3 * sol.sz.sz );

    [ Vm ] = edgedetect_FDxy( phi, sFDxy );

    %===================

%     lam.x = lam.x + 0.9999 * ( Vs.x - Vm.x );
%     lam.y = lam.y + 0.9999 * ( Vs.y - Vm.y );
    lam.x = lam.x + 1.0 * ( Vs.x - Vm.x );
    lam.y = lam.y + 1.0 * ( Vs.y - Vm.y );
