%

%{

clear; close all; 
restoredefaultpath
code_path = '~/Documents/Science/Matlab/Code/cdi';
cd( code_path );
addpath( genpath( code_path ));

clear; close all; testing_sparsedictlearn_phaseret

%}


%==================================================================================================

% new random seed
rng shuffle;

%==================================================================================================

% load data

rsdata = '~/Documents/Science/Matlab/Code/cdi/run/misctesting2DTPA/simPRexpt_2DTPA.mat';
load( rsdata );

expt.paths.rsdata = rsdata;
clear rsdata;

% sol.it.exwv = 1;

%==================================================================================================

sol.plot.xaxis = expt.csys.z2.dc * (1 : sol.sz.c) * 1e6;
sol.plot.xaxis = sol.plot.xaxis - mean( sol.plot.xaxis );
sol.plot.yaxis = expt.csys.z2.dr * (1 : sol.sz.r) * 1e6; 
sol.plot.yaxis = sol.plot.yaxis - mean( sol.plot.yaxis );

%==================================================================================================

tmp0 = make_gaussian( sol.sz.sz, [ 0.5 * sol.sz.r + 1, 0.5 * sol.sz.c + 1 ], [ 1.999 * sol.sz.r, 1.999 * sol.sz.c ]);
sol.measLPF = 0 + 1 * fftshift( tmp0 );

% sol.probe.S = 0 + 1 * make_rectangle( sol.sz.sz, [ 0.8 * sol.sz.r, 0.8 * sol.sz.c ]);
sol.probe.S = 0 + 1 * make_ellipsoid( sol.sz.sz, [ 0.9 * sol.sz.c, 0.9 * sol.sz.r ]);
[ sol.probe.S ] = lpf_gauss( sol.probe.S, 440.01 * sol.sz.sz );

clear( 'tmp0' )

%==================================================================================================


if 0 || ( sol.it.exwv == 1 )
    
    sol.ss = 1;
    
    [ expt.phiT, ~ ] = enforce_2DTPAsposview( expt.probe.P, expt.sample.TF, expt.sample.vs, expt.spos.rs( sol.ss, : ), 'subpx' );
    
    [ sol.phiB, ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, expt.sample.vs, expt.spos.rs( sol.ss, : ), 'subpx' );
%     sol.phiB        = enforce_2DTPAmeas( sol.phiB, expt.meas.SI( sol.ss ), sol.measLPF, sol );

%     sol.synthesisDL = setup_synthesis_dictionarylearning( sol.sz.sz, expt );
%     sol.analysisDL  = setup_analysis_dictionarylearning( sol.sz.sz );
%     sol.sparseFDxy  = setup_edgedetect_forwarddiff2Dxy( [], sol.sz );   
    
%     [U,S,V] = svd( Ip, 'econ' );
%     [U,S,V] = svd( sol.analysisDL.D );
%     Sinv = 1 ./ S; 
%     Sinv( Sinv > 1e6 ) = 0;
%     Sinv2 = pinv(S);
%     D = Ca * ( V * pinv(S) * U' );
    

%     sol.probe.P = expt.probe.P;

end

% sol.phiB = expt.phiT;

% sol.phiB = abs( sol.phiB );

sol.spos.updateorder = sol.ss;

%==================================================================================================

% CENTER THE EXIT WAVES:


% figure; imagesc( abs( sol.phiB ))
% 
% [ com ] = centerOfMass( abs( sol.phiB ));
% 
% sol.phiB = circshift( sol.phiB, -1 * round( com - 0.5 * sol.sz.sz - 1 ));
% 
% figure; imagesc( abs( sol.phiB ))


% [~, Ic ] = max( sum( abs( sol.phiB ), 1 ));
% [~, Ir ] = max( sum( abs( sol.phiB ), 2 ));
% 
% sol.phiB = circshift( sol.phiB, -1 * round( [ Ir, Ic ] - 0.5 * sol.sz.sz - 1 ));
% 
% figure; imagesc( abs( sol.phiB ))









%==================================================================================================

[ sol.sDL_overcomplete ] = setup_dictSYNTHESIS_overcomplete( sol.sz.sz, expt );

sol.sDL_overcomplete.mean0 = logical( 0 ); %#ok<LOGL>

% number of nonzeros to keep in the sparse synthesis code matrix columns
sol.sDL_overcomplete.sparse_lvl = 0.05;

sol.sDL_overcomplete.Cnnz = round( sol.sDL_overcomplete.Na * sol.sDL_overcomplete.sparse_lvl );              
if sol.sDL_overcomplete.Cnnz <= 0, sol.sDL_overcomplete.Cnnz = 1; end

% sol.sDL_overcomplete.Cnnz = 1;

%=======================

[ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_overcomplete, sol.sz.sz );
sol.sDL_overcomplete.C = sDL_sparsecode_update( IpT, sol.sDL_overcomplete );

%==================================================================================================

[ sol.sDL_orthogonal ] = setup_dictSYNTHESIS_orthogonal( sol.sz.sz, expt );

sol.sDL_orthogonal.mean0 = logical( 0 ); %#ok<LOGL>

% number of nonzeros to keep in the sparse synthesis code matrix columns
sol.sDL_orthogonal.sparse_lvl = 0.05;

sol.sDL_orthogonal.Cnnz = 0 + 1 * round( sol.sDL_orthogonal.Na * sol.sDL_orthogonal.sparse_lvl );              
if sol.sDL_orthogonal.Cnnz <= 0, sol.sDL_orthogonal.Cnnz = 1; end

[ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal );

%=======================

sol.sDL_orthogonal.ADMMlam = complex( zeros( sol.sz.sz, 'single' ));

sol.sDL_orthogonal.thresh = 'soft';
sol.sDL_orthogonal.thresh = 'hard';

%=======================

% sol.sDL_orthoRE = sol.sDL_orthogonal;
% sol.sDL_orthoIM = sol.sDL_orthogonal;
% 
% [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
% sol.sDL_orthoRE.C = sDL_sparsecode_update( real( IpT ), sol.sDL_orthoRE );
% sol.sDL_orthoIM.C = sDL_sparsecode_update( imag( IpT ), sol.sDL_orthoIM );

%==================================================================================================


















% % sDLp = sol.sDL_orthogonal;
% sDLp = sol.sDL_overcomplete;
% 
% 
% 
% % get true exit wave image patch training sets
% [ IpT, ~ ] = image2trainingsetpatches( expt.phiT, sDLp, sol.sz.sz );
% 
% [ f2, f1, W ] = trainingsetpatches2image( IpT, sDLp, sol.sz.sz );
% 
% 
% % re_IpT = real( IpT );
% % im_IpT = imag( IpT );  
% re_IpT = abs( IpT );
% im_IpT = angle( IpT );  
% 
% ll = 1;
% 
% 
% figure; 
% imagesc( abs( sDLp.D' * sDLp.D )); colormap gray; colorbar
% figure; 
% plot_dictionnary( sDLp, [ 10, 10 ] );
% 
% for kk = 1 : 20
%     
%     kk
%     
% %     %==================
% %     
% %     if ( mod( kk, 1 ) == 0 )
% %         
% %         tmp1 = sDLp.D * sDLp.C;
% %         [ tmp0 ] = trainingsetpatches2image( tmp1, sDLp, sol.sz.sz );
% %         [ Ip, ~ ] = image2trainingsetpatches( expt.phiT, sDLp, sol.sz.sz );
% %     
% %         figure; 
% %         subplot(121); imagesc( abs( tmp0 )); colormap gray; axis square
% %         subplot(122); imagescHSV( tmp0 ); axis square
% % 
% %     end
% %     
% %     %==================
% 
%             
%     
%     
%     
%     
% 
%     %==================
%     
%     if ( mod( kk, 1 ) == 0 ), 
%         
%         sDLp.D = sDL_dictionary_update( IpT, sDLp ); 
%         
% %         sDLp.D = real( sDLp.D );
% figure; 
% plot_dictionnary( sDLp, [ 10, 10 ] );
%     5;
%     
%     end
%     sDLp.C = sDL_sparsecode_update( IpT, sDLp );
%     
%     
% % sol.sDL_orthoRE.C = sDL_sparsecode_update( re_IpT, sol.sDL_orthoRE );
% % sol.sDL_orthoIM.C = sDL_sparsecode_update( im_IpT, sol.sDL_orthoIM );
% % sDLp.C =  sol.sDL_orthoRE.C + 1i *  sol.sDL_orthoIM.C;
%     
%     
%     
%     
%     
% %             sDLp.D = sDL_dictionary_update( IpT, sDLp );
% %             sol.sDL_orthoRE.D = abs( sDLp.D );
% %             sol.sDL_orthoIM.D = angle( sDLp.D );   
% 
% %             sDLp.D = sDL_dictionary_update( IpT, sDLp );
% %             sol.sDL_orthoRE.D = abs( sDLp.D );
% %             sol.sDL_orthoIM.D = angle( sDLp.D );    
% %     
% 
% %             sDLp.D = sDL_dictionary_update( IpT, sDLp );
% %             sol.sDL_orthoRE.D = real( sDLp.D );
% %             sol.sDL_orthoIM.D = imag( sDLp.D );    
% 
% %             sol.sDL_orthoRE.D = sDL_dictionary_update( re_IpT, sol.sDL_orthoRE );
% %             sol.sDL_orthoIM.D = sDL_dictionary_update( im_IpT, sol.sDL_orthoIM );
% %             sDLp.D = sol.sDL_orthoRE.D .* exp( 1i * sol.sDL_orthoIM.D ); 
% 
% %             sol.sDL_orthoRE.D = sDL_dictionary_update( re_IpT, sol.sDL_orthoRE );
% %             sol.sDL_orthoIM.D = sDL_dictionary_update( im_IpT, sol.sDL_orthoIM );
% %             sDLp.D = sol.sDL_orthoRE.D + 1i * sol.sDL_orthoIM.D; 
%             
% %             figure; 
% %             plot_dictionnary( sDLp, [ sDLp.wy, sDLp.wx ] );
% %             close all;
%             
% %             sDLp.C = sDL_sparsecode_update( IpT, sDLp );
% 
% %             sDLp.C = sDL_sparsecode_update( IpT, sDLp );
% %             sol.sDL_orthoRE.C = real( sDLp.C );
% %             sol.sDL_orthoIM.C = imag( sDLp.C );
%             
% %             sol.sDL_orthoRE.C = sDL_sparsecode_update( re_IpT, sol.sDL_orthoRE );
% %             sol.sDL_orthoIM.C = sDL_sparsecode_update( im_IpT, sol.sDL_orthoIM );
% %             sDLp.C = sol.sDL_orthoRE.C + 1i * sol.sDL_orthoIM.C; 
% 
% %             sol.sDL_orthoRE.C = sDL_sparsecode_update( re_IpT, sol.sDL_orthoRE );
% %             sol.sDL_orthoIM.C = sDL_sparsecode_update( im_IpT, sol.sDL_orthoIM );
% %             sDLp.C = sol.sDL_orthoRE.C .* exp( 1i * sol.sDL_orthoIM.C ); 
%             
%     %==================
%     
%     if ( mod( kk, 5 ) == 0 )
%         
%         tmp1 = sDLp.D * sDLp.C;
%         [ tmp0 ] = trainingsetpatches2image( tmp1, sDLp, sol.sz.sz );
%         [ Ip, ~ ] = image2trainingsetpatches( expt.phiT, sDLp, sol.sz.sz );
%         
%         if ( mod( kk, 50 ) == 0 )
%             
%             figure; 
%             subplot(121); imagesc( abs( tmp0 )); colormap gray; axis square
%             subplot(122); imagescHSV( tmp0 ); axis square
%         
%         end
%         
%         err1( ll ) = norm( Ip - tmp1, 'fro' );
%         err2( ll ) = norm( expt.phiT - tmp0, 'fro' );
%         ll = ll + 1;
%     
%     end
% 
% 
%     
%     %==================
%     
% 
% 
% end
% 
% tmp1 = sDLp.D * sDLp.C;
% [ tmp2 ] = trainingsetpatches2image( tmp1, sDLp, sol.sz.sz );
% 
% figure; 
% subplot(121); imagesc( abs( tmp2 )); colormap gray; axis square
% subplot(122); imagescHSV( tmp2 ); axis square
% 
% figure; 
% subplot(121); imagesc( abs( expt.phiT )); colormap gray; axis square
% subplot(122); imagescHSV( expt.phiT ); axis square
% 
% figure; 
% subplot(211); semilogy( err1, '-o' ); title('|| Ip - Ds Cs ||_F')
% subplot(212); semilogy( err2, '-o' ); title('|| \phi_T - \phi^{(k)} ||_F')
% % subplot(313); semilogy( err3, '-o' ); title('|| Cs ||_0')
% 
% figure; 
% imagesc( abs( sDLp.D' * sDLp.D )); colormap gray; colorbar
% figure; 
% plot_dictionnary( sDLp, [ 10, 10 ] );
% 
% figure; 
% imagesc( abs( sDLp.C )~=0 )
% colormap gray
% axis square
% 
% 
% norm( tmp2 - expt.phiT, 'fro' )
% 
% 
% % close all;
% return












%==================================================================================================

sol.sparseFDxy  = setup_edgedetect_forwarddiff2Dxy( [], sol.sz );  

% sol.sparseFDxy.sparselvl = round( sol.sz.rc * 0.12 ); 
% 
% % sol.sparseFDxy.nzmask = 'anisotropic';
% sol.sparseFDxy.nzmask = 'isotropic';
% % sol.sparseFDxy.thresh = 'soft';
% sol.sparseFDxy.thresh = 'hard';


% lagrangian multipliers for sparse forward differences
lam.x = 0 * sol.phiB;
lam.y = lam.x;
[ Vm ] = edgedetect_FDxy( sol.phiB, sol.sparseFDxy );


[ Vt ] = edgedetect_FDxy( expt.phiT, sol.sparseFDxy );
[ V ] = iedgedetect_FDxy( Vt, sol.sparseFDxy, sol.sz );




% V_fft = ( sol.sparseFDxy.Dx_fft_conj .* fft2( Vt.x ) ) / sol.sz.sqrt_rc;
% V_fft = V_fft .* abs( sol.sparseFDxy.Dx_fft_conj ).^2;

V_fft = fft2( Vt.x ) / sol.sz.sqrt_rc;
V_fft = V_fft ./ ( 1e-7 + conj( sol.sparseFDxy.Dx_fft_conj ));

V_fft( 1 ,1 ) = Vt.q0 / sol.sz.sqrt_rc; 

V2 = ifft2( V_fft ) * sol.sz.sqrt_rc;

figure; imagesc( abs(expt.phiT) ); daspect([1 1 1]); colormap jet
figure; imagesc( abs(V) ); daspect([1 1 1]); colormap jet
figure; imagesc( abs(V2) ); daspect([1 1 1]); colormap jet





%==================================================================================================

[ sol.aDL ] = setup_dictANALYSIS_overcomplete( sol.sz.sz, expt );


% number of nonzeros to keep in the sparse synthesis code matrix columns
sol.aDL.sparse_lvl = 0.09;
sol.aDL.Cnnz = 0 + 1 * round( sol.aDL.Na * sol.aDL.sparse_lvl );              
if sol.aDL.Cnnz <= 0, sol.aDL.Cnnz = 1; end

[ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.aDL, sol.sz.sz );
sol.aDL.C = aDL_sparsecode_update( IpT, sol.aDL );

aDLlam = 0 * sol.aDL.C;







%==================================================================================================

% sol.synthesisDL.C = eye( sol.synthesisDL.Na, sol.synthesisDL.Np, 'single' ); 

% sol.phiB = expt.phiT;
% [ IpT, ~ ] = image2trainingsetpatches( expt.phiT, sol.synthesisDL, sol.sz.sz );
% sol.synthesisDL.C = sDL_sparsecode_update( IpT, sol.synthesisDL );
% sol.synthesisDL.D = sDL_dictionary_update( IpT, sol.synthesisDL.C, sol.synthesisDL.D );

% [ sol.phiB ] = trainingsetpatches2image( sol.synthesisDL.D * sol.synthesisDL.C, sol.synthesisDL, sol.sz.sz );
% [ sol.phiB, sol.synthesisDL ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );


Vs = Vm;

for ii = 1 : 4000
    
    fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
    if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end
    
    %==============================================================================================
    %------ PHASE RETRIEVAL ON EXIT WAVE USING FORWARD DIFFERENCE SPARSIFYING TRANSFORMATION ------
    %==============================================================================================

    
%     [ sol.phiB, Vm, Vs, lam ] = ADMMupdate_exwv_sparseFDxy_meas( Vm, lam, sol.sparseFDxy, expt.meas.SI( sol.ss ), sol );
%     [ sol.phiB ] = RAARupdate_exwv_sparseFDxy_meas( sol.phiB, sol.sparseFDxy, expt.meas.SI( sol.ss ), sol );
%     [ Vs, sol.phiB ] = RAARupdate_exwv_sparseFDxy_meas_V2( Vs, sol.sparseFDxy, expt.meas.SI( sol.ss ), sol );



[ sol.phiB, sol.aDL.C, aDLlam, aDL ] = ADMMupdate_exwv_aDL_C_meas( sol.aDL.C, aDLlam, aDL, expt.meas.SI( sol.ss ), sol );


    %==============================================================================================
    %----------- PHASE RETRIEVAL ON EXIT WAVE USING ANALYSIS SPARSE DICTIONARY LEARNING -----------
    %==============================================================================================

%     [ sol.phiB, VmB, VsB, lamB, sol.analysisDL ] = ADMMupdate_exwv_aDL_meas( VmB, VsB, lamB, sol.analysisDL, expt.meas.SI( sol.ss ), sol );
    


    % TRY HIO, ER, OR RAAR ON EXITWAVE
    
    
    
    
    
    if ( mod( sol.it.exwv, 1e99 ) == 0 )

         

%         sol.phiB = enforce_2DTPAmeas( sol.phiB, expt.meas.SI( sol.ss ), sol.measLPF, sol );
%         [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
%         sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal );
% %         sol.sDL_orthogonal.D = sDL_dictionary_update( IpT, sol.sDL_orthogonal.C, sol.sDL_orthogonal.D );
%         [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz.sz );




% %         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         [ sol.phiB, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
% %         
%         [ sol.phiB, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
% %         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );


%         [ phiC, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );  
% %         [ phiD, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         [ phiD, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         sol.phiB = 0.5 * ( phiC + phiD );
        


        
%         [ phiC, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         [ phiD, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         sol.phiB = 0.5 * ( phiC + phiD );
        
        
%         
%         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );





% [ sol.phiB, Vs, sol.ADMMlam, sol.sDL_orthogonal ] = ADMMupdate_exwv_DsCs_meas( sol.phiB, sol.ADMMlam, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
% 





% [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_DsCs_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );


[ sol.phiB, sol.sDL_orthogonal.C, ~ ] = RAARupdate_exwvmeas_sDLsparsecode( sol.phiB, sol.sDL_orthogonal.C, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );

sol.phiB = abs( sol.phiB );

if (mod( sol.it.exwv, 1 ) == 0 )
    
    for uu = 1 : 1
        
        [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
        sol.sDL_orthogonal.D = sDL_dictionary_update( IpT, sol.sDL_orthogonal );
%         sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal );
%         [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz.sz );
%         sol.phiB = enforce_2DTPAmeas( sol.phiB, expt.meas.SI( sol.ss ), sol.measLPF, sol );
    
    end
    
end




















% %         [ tmp0 ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz.sz );
% %         sol.phiB = enforce_2DTPAmeas( tmp0, expt.meas.SI( sol.ss ), sol.measLPF, sol );
% 
%         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_DsCs_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
%         
%         sol.phiB = abs( sol.phiB );
%         
%         [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
%         
%         if (mod( sol.it.exwv, 1 ) == 0 )
%             
% %             [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
%             sol.sDL_orthogonal.D = sDL_dictionary_update( IpT, sol.sDL_orthogonal );
%             sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal );
% %             [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz.sz );
%             
% % sol.sDL_orthogonal.D = abs( sol.sDL_orthogonal.D );
% % sol.sDL_orthogonal.C = abs( sol.sDL_orthogonal.C );
% 
%         end

        
        
        
        
        
        
        
        
        
        
        
        
% %         [ tmp0 ] = trainingsetpatches2image( sol.sDL_orthoRE.D * sol.sDL_orthoRE.C, sol.sDL_orthoRE, sol.sz.sz );
% %         sol.phiB = enforce_2DTPAmeas( tmp0, expt.meas.SI( sol.ss ), sol.measLPF, sol );
% 
%         [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_DsCs_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
% %         [ sol.phiB, sol.sDL_orthoRE, sol.sDL_orthoIM ] = RAARupdate_exwv_ReImDsCs_meas( sol.phiB, sol.sDL_orthoRE, sol.sDL_orthoIM, expt.meas.SI( sol.ss ), sol );
%         
%         [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
%         
%         if (mod( sol.it.exwv, 1 ) == 0 )
%             
%             re_IpT = real( IpT );
%             im_IpT = imag( IpT );   
%             
%             sol.sDL_orthoRE.D = sDL_dictionary_update( re_IpT, sol.sDL_orthoRE );
%             sol.sDL_orthoIM.D = sDL_dictionary_update( im_IpT, sol.sDL_orthoIM );
%             sol.sDL_orthogonal.D = sol.sDL_orthoRE.D + 1i * sol.sDL_orthoIM.D; 
% %             
% %             figure; 
% %             plot_dictionnary( sol.sDL_orthogonal, [ sol.sDL_orthogonal.wy, sol.sDL_orthogonal.wx ] );
% %             close all;
% 
%             sol.sDL_orthoRE.C = sDL_sparsecode_update( re_IpT, sol.sDL_orthoRE );
%             sol.sDL_orthoIM.C = sDL_sparsecode_update( im_IpT, sol.sDL_orthoIM );
%             sol.sDL_orthogonal.C = sol.sDL_orthoRE.C + 1i * sol.sDL_orthoIM.C; 
%             
% 
% 
% 
% %             [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthoRE.D * sol.sDL_orthoRE.C, sol.sDL_orthoRE, sol.sz.sz );
%         
%         end
%         

        
        
        
        
        
        
        
        
        
        
        
        
        
% %         [ sol.phiB, sol.sDL_overcomplete ] = RAARupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_overcomplete, expt.meas.SI( sol.ss ), sol );
% % %         [ sol.phiB, sol.sDL_overcomplete ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_overcomplete, expt.meas.SI( sol.ss ), sol );
% 
% 
% 
%         [ sol.phiB, sol.sDL_overcomplete ] = RAARupdate_exwv_DsCs_meas( sol.phiB, sol.sDL_overcomplete, expt.meas.SI( sol.ss ), sol );
%         
%         [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_overcomplete, sol.sz.sz );
%         
%         if (mod( sol.it.exwv, 1 ) == 0 )
%             
% %             [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_overcomplete, sol.sz.sz );
%             sol.sDL_overcomplete.D = sDL_dictionary_update( IpT, sol.sDL_overcomplete );
%             sol.sDL_overcomplete.C = sDL_sparsecode_update( IpT, sol.sDL_overcomplete );
% %             [ sol.phiB ] = trainingsetpatches2image( sol.sDL_overcomplete.D * sol.sDL_overcomplete.C, sol.sDL_overcomplete, sol.sz.sz );
%         
%         end












%         sol.phiB = abs( sol.phiB );
%         sol.sDL_orthogonal.D = abs( sol.sDL_orthogonal.D );
%         sol.sDL_orthogonal.C = abs( sol.sDL_orthogonal.C );



        if (mod( sol.it.exwv, 10e9 ) == 0 )
            
%             [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
% %             sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal );
%             sol.sDL_orthogonal.D = sDL_dictionary_update( IpT, sol.sDL_orthogonal.C, sol.sDL_orthogonal.D );
%             [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz.sz );
            
%             [ sol.phiB, sol.sDL_orthogonal ] = RAARupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );   
    
            for jj = 1 : 5
                [ sol.phiB, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
                [ sol.phiB, sol.sDL_orthogonal ] = ERupdate_exwv_sDL_D_meas( sol.phiB, sol.sDL_orthogonal, expt.meas.SI( sol.ss ), sol );
            end
                
        end
        
%         [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.sDL_orthogonal, sol.sz.sz );
%         sol.sDL_orthogonal.C = sDL_sparsecode_update( IpT, sol.sDL_orthogonal.C, sol.sDL_orthogonal.D, sol.sDL_orthogonal.Cnnz, sol.sDL_orthogonal.thresh );
%         sol.phiB = sol.phiB .* sol.probe.S;

%         [ sol.phiB ] = trainingsetpatches2image( sol.sDL_orthogonal.D * sol.sDL_orthogonal.C, sol.sDL_orthogonal, sol.sz );


    end
    
    %==============================================================================================
    
%     if (mod( sol.it.exwv, 1e9 ) == 0 )
% 
% %         [ sol.phiB, sol.analysisDL ] = RAARupdate_exwv_aDL_C_meas( sol.phiB, sol.analysisDL, expt.meas.SI( sol.ss ), sol );
%         [ sol.phiB, sol.analysisDL ] = ERupdate_exwv_aDL_C_meas( sol.phiB, sol.analysisDL, expt.meas.SI( sol.ss ), sol );
%         
% %         [ sol.phiB, sol.analysisDL ] = RAARupdate_exwv_aDL_D_meas( sol.phiB, sol.analysisDL, expt.meas.SI( sol.ss ), sol );
%         [ sol.phiB, sol.analysisDL ] = ERupdate_exwv_aDL_D_meas( sol.phiB, sol.analysisDL, expt.meas.SI( sol.ss ), sol );
% 
% 
%     end
    
    %==============================================================================================
    
    
    
    
    
    
    %==============================================================================================
    
    if mod( sol.it.exwv, 200 ) == 0 

        sol.meas_error( sol.it.mtot ) = norm( expt.meas.SI( sol.ss ).D - abs( fft2( fftshift( sol.phiB ))) / sol.sz.sqrt_rc  , 'fro' );
        sol.it.metr( sol.it.mtot ) = sol.it.exwv;
        sol.it.mtot = sol.it.mtot + 1;

        figure( 666 ); 
        set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
%         subaxis( 2, 3, 1, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 ); 
        subplot( 2, 3, [ 1, 4 ] )
        plot( sol.it.metr, sol.meas_error, 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] ); 
%         semilogy( sol.it.metr, sol.meas_error, 'linewidth', 2, 'color', [ 0.0, 0.8, 0.0 ] ); 
        
%         subaxis( 2, 3, 2, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
        subplot( 2, 3, 2 )
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( sol.phiB ))); 
        colorbar
        xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
        title( [ 'abs exwv', num2str( sol.it.exwv, ', it = %d' )])
        
%         subaxis( 2, 3, 3, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
        subplot( 2, 3, 3 )
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( expt.phiT ))); 
        colorbar
        xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
        title( 'TRUE abs exwv' )
                
        

%         subaxis( 2, 3, 5, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
        subplot( 2, 3, 5 )
        imagescHSV( log10( 1 + abs( sol.phiB )) .* exp( 1i * angle( sol.phiB )), sol.plot ); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
        title( 'HSV abs exwv' )
        
%         subaxis( 2, 3, 6, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
        subplot( 2, 3, 6 )
        imagescHSV( log10( 1 + abs( expt.phiT )) .* exp( 1i * angle( expt.phiT )), sol.plot ); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
        title( 'HSV TRUE  abs exwv' )
        
        pause( 0.02 )
        
        export_fig( num2str( [ sol.it.exwv, sol.ss ], 'exitwave_%d_spos-%d.jpg' ), '-r90.0' )
        close all;
        
    end
    
    %==================================
    % EXITWAVE ITERATION COUNTER
    %==================================
    
    sol.it.exwv = sol.it.exwv + 1;
    
end




clearvars -except expt sol

% sol.phi = []; sol = rmfield( sol, 'phi' );
% sol.phiOLD = []; sol = rmfield( sol, 'phiOLD' );
% sol.unmeas = []; sol = rmfield( sol, 'unmeas' );

% save data
fprintf('\n========================================================================================================'); 
fprintf('\nSAVING Simulated Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed...'); 
save( expt.paths.rsdata, '*' );
fprintf('done saving data!\n'); 
fprintf('========================================================================================================\n\n'); 







%==================================================================================================



% N = 50; 
% 
% A = rand( N, N ) + 1i * rand( N, N ); 
% 
% norm( A, 'fro' )^2
% 
% 
% As = 0;
% 
% for ii = 1 : size( A, 2 )
%     
%     As = As + norm( A( :, ii ),2)^2;
%     
% end
% 
% As




%{


r = 0.5;

R1 = [ [ 1, 0, 0, 0, 0, 0, 0, 0, 0 ]; [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ]; [ 0, 0, 0, 1, 0, 0, 0, 0, 0 ];  [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ]; ];

R1' * R1
inv( R1' * R1 + r * eye(9))



R2 = [ [ 0, 0, 0, 1, 0, 0, 0, 0, 0 ]; [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ]; [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ];  [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ]; ];

R2' * R2

inv( R2' * R2 +  1.0 * eye(9))

inv( R1' * R1 + R2' * R2 +  r * eye(9))



R3 = [ [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ]; [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ]; [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ];  [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ]; ];

R3' * R3

inv( R3' * R3 +  1.0 * eye(9))



inv( R1' * R1 + R2' * R2 + R3' * R3 +  r * eye(9))






clear

R = [ [ 1, 0, 0, 0, 0 ];
      [ 0, 1, 0, 0, 0 ];
      [ 0, 0, 0, 1, 0 ];
      [ 0, 0, 0, 0, 1 ]; ];  

R1 =[ zeros( 4, 0 ), R, zeros( 4, 7 ) ];

R2 = [ zeros( 4, 3 ), R, zeros( 4, 4 ) ];

R3 = [ zeros( 4, 6 ), R, zeros( 4, 1 ) ];

R4 = [ zeros( 4, 1 ), R, zeros( 4, 6 ) ];

R5 = [ zeros( 4, 4 ), R, zeros( 4, 3 ) ];

R6 = [ zeros( 4, 7 ), R, zeros( 4, 0 ) ];

R1' * R1 + R2' * R2 + R3' * R3 + R4' * R4 + R5' * R5 + R6' * R6






R1b = [ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ];  
        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ]; ];

norm( R1 - R1b, 'fro')


R2b = [ [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ];  
        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]; ];

norm( R2 - R2b, 'fro')

R3b = [ [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ];  
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ]; ];

norm( R3 - R3b, 'fro')

R4b = [ [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ];  
        [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ]; ];

norm( R4 - R4b, 'fro')


%}





%{ 

D0 = sol.synthesisDL;

% subtract off the mean of the image patches


% % get true exit wave image patch training set
% [ IpT, ~ ] = image2trainingsetpatches( expt.phiT, sol.synthesisDL, sol.sz.sz );
% 
% % compute initial sparse code using initial image patch training set
% sol.synthesisDL.C = sDL_sparsecode_update( IpT, sol.synthesisDL.C, sol.synthesisDL.D, sol.synthesisDL.Cnnz, sol.synthesisDL.thresh );
% % sol.synthesisDL.D = sDL_dictionary_update( IpT, sol.synthesisDL.C, sol.synthesisDL.D );
% 
% % tmp1 = sol.synthesisDL.D * sol.synthesisDL.C;
% % [ tmp0 ] = trainingsetpatches2image( tmp1, sol.synthesisDL, sol.sz.sz );
% % figure; imagesc( abs( tmp0 ))

a = 1.0;
sol.phiB = a * sol.phiB + ( 1 - a ) * expt.phiT;

for ii = 1 : 1000
    
    ii
    
    %===================
    
    tmp1 = sol.synthesisDL.D * sol.synthesisDL.C;
    [ tmp0 ] = trainingsetpatches2image( tmp1, sol.synthesisDL, sol.sz.sz );
    [ tmp2, ~ ] = image2trainingsetpatches( expt.phiT, sol.synthesisDL, sol.sz.sz );
    
    err1( ii ) = norm( tmp2 - tmp1, 'fro' ) + 0.0 * sum( sol.synthesisDL.C( : ) ~= 0 );
    err2( ii ) = norm( expt.phiT - tmp0, 'fro' );
    
    clear( 'tmp0', 'tmp1', 'Ip' )
    
    %===================
    
%     [ sol.phiB, sol.synthesisDL ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );
    [ sol.phiB, sol.synthesisDL ] = RAARupdate_exwv_sDL_C_meas( sol.phiB, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );
    
%     [ ~, sol.synthesisDL ] = ERupdate_exwv_sDL_D_meas( expt.phiT, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );
%     [ ~, sol.synthesisDL ] = RAARupdate_exwv_sDL_D_meas( expt.phiT, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );
    
    %===================
    
    sol.phiB = sol.phiB .* sol.probe.S;
    
%     sol.phiB = sol.phiB * 15 / max( abs( sol.phiB( : )));
    
%     max_abs = 15;
%     tmp4 = abs( sol.phiB ) > max_abs;
%     sol.phiB( tmp4 ) = max_abs * exp( 1i * angle( sol.phiB( tmp4 )));
    
%     sol.phiB = lpf_gauss( sol.phiB, 0.1 * sol.sz.sz );

    %===================
    
    if mod( ii, 10 ) == 0 
        
        [ sol.phiB, sol.synthesisDL ] = ERupdate_exwv_sDL_C_meas( sol.phiB, sol.synthesisDL, expt.meas.SI( sol.ss ), sol );
        [ IpT, ~ ] = image2trainingsetpatches( sol.phiB, sol.synthesisDL, sol.sz.sz );
        sol.synthesisDL.C = sDL_sparsecode_update( IpT, sol.synthesisDL.C, sol.synthesisDL.D, sol.synthesisDL.Cnnz, sol.synthesisDL.thresh );
%         sol.synthesisDL.D = sDL_dictionary_update( IpT, sol.synthesisDL.C, sol.synthesisDL.D );

    end
    
    %===================
    
%     figure; imagesc( abs( sol.synthesisDL.C ~= 0 ))
%     
%     colsD_rm = ~any( sol.synthesisDL.D, 1 );
%     sum( colsD_rm )
    
%     tmp1 = sol.synthesisDL.D * sol.synthesisDL.C;
%     [ tmp0 ] = trainingsetpatches2image( tmp1, sol.synthesisDL, sol.sz );
%     [ Ip, ~ ] = image2trainingsetpatches( expt.phiT, sol.synthesisDL, sol.sz );
%     
%     err1( ii ) = norm( Ip - tmp1, 'fro' ) + 0.1 * sum( sol.synthesisDL.C( : ) ~= 0 );
%     err2( ii ) = norm( expt.phiT - tmp0, 'fro' );
    
    %===================
    
    
    if mod( ii, 10 ) == 0 

        figure( 666 ); 
        set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        subaxis( 2, 3, 1, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 ); 
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( sol.phiA ))); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
        
        subaxis( 2, 3, 2, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( sol.phiB ))); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
        
        subaxis( 2, 3, 3, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( expt.phiT ))); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
          
        

        subaxis( 2, 3, 4, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
        imagescHSV( log10( 1 + abs( sol.phiA )) .* exp( 1i * angle( sol.phiA )), sol.plot ); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])

        subaxis( 2, 3, 5, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
        imagescHSV( log10( 1 + abs( sol.phiB )) .* exp( 1i * angle( sol.phiB )), sol.plot ); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
        title( num2str( ii, '%d' ))
        
        subaxis( 2, 3, 6, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
        imagescHSV( log10( 1 + abs( expt.phiT )) .* exp( 1i * angle( expt.phiT )), sol.plot ); 
        xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
        
        pause( 0.02 )
        
        export_fig( num2str( [ ii, sol.ss ], 'exitwave_%d_spos-%d.jpg' ), '-r90.0' )
        close all;
        
    end
    
    
    
    
    
end

tmp1 = sol.synthesisDL.D * sol.synthesisDL.C;
[ tmp0 ] = trainingsetpatches2image( tmp1, sol.synthesisDL, sol.sz.sz );

figure; 
subplot(221); imagescHSV(expt.phiT); daspect([1 1 1])
subplot(222); imagescHSV(tmp0); daspect([1 1 1])
subplot(223); imagesc(abs(expt.phiT)); daspect([1 1 1]); colorbar
subplot(224); imagesc(abs(tmp0)); daspect([1 1 1]); colorbar
% subplot(224); imagesc(abs(tmp0 - expt.phiT)); daspect([1 1 1]); colorbar

figure; semilogy( err1 ); title('|| Ip - Ds Dc ||_F')
figure; semilogy( err2 ); title('|| \phi_T - \phi^{(k)} ||_F')


Df = sol.synthesisDL;
colsD_rm = ~any( sol.synthesisDL.D, 1 );

% figure; imagesc( sol.synthesisDL.D == 0 )

Df.D( :, colsD_rm ) = []; 
D0.D( :, colsD_rm ) = [];

figure; plot_dictionnary( D0, [ 4, 4 ] );
figure; plot_dictionnary( Df, [ 4, 4 ] );

% close all;
return

%}




% % [ Q1, R1 ] = mgsog( sol.synthesisDL.D );
% % [ Q2, R2 ] = gsog( sol.synthesisDL.D );
% % [ Q3, R3 ] = gson( sol.synthesisDL.D );
% % [ Q4, R4 ] = mgsog( sol.synthesisDL.D );
% % [ Q5, R5 ] = mgsog( sol.synthesisDL.D );
% % [ Q6, R6 ] = qr( sol.synthesisDL.D );
% % % Ds = orth( sol.synthesisDL.D );
% 
% 
% 
% [ U, S, V ] = svd( sol.synthesisDL.D, 0 );
% % S( S ~= 0 ) = 1;  % S( diag ) = 1;
% S = eye( sol.synthesisDL.Ni, sol.synthesisDL.Na );
% Ds = U * S * V';
% 
% Vt = V';
% Vt = Vt( 1 : sol.synthesisDL.Ni, : ); 
% Ds2 = U * Vt;
% 
% figure; imagesc( abs( Ds' * Ds )); colorbar
% figure; imagesc( abs( Ds2' * Ds2 )); colorbar
% figure; imagesc( abs( pinv(Ds) * Ds )); colorbar
% 
% 
% [ U, E ] = eig( sol.synthesisDL.D' * sol.synthesisDL.D ); 
% % Ds3 = sol.synthesisDL.D * U ./ ( 1e-7 + repmat( sqrt( diag( E )).' , [ sol.synthesisDL.Ni, 1 ] ) ); 
% Ds3 = sol.synthesisDL.D * U; 
% 
% figure; imagesc( log10( 1 + 10^5 * abs( Ds3' * Ds3 ))); colorbar
% figure; imagesc( log10( 1 + 10^5 * abs( E ))); colorbar
% 
% 
% [ U, E ] = eig( sol.synthesisDL.D * sol.synthesisDL.D' ); 
% Ds4 = U * sol.synthesisDL.D; 
% 
% % figure; imagesc( log10( 1 + 10^0 * abs( Ds4' * Ds4 ))); colorbar
% figure; imagesc( abs( Ds4' * Ds4 )); colorbar
% 
% 
% [ U, S, V ] = svd( sol.synthesisDL.D, 0 );
% S( S ~= 0 ) = 1;
% synthesisDL.D = U * S * V';
% 
% % [ U1, E1 ] = eig( sol.synthesisDL.D' * sol.synthesisDL.D ); 
% % tmp0 = sol.synthesisDL.D * U1;
% % 
% % [ U2, E2 ] = eig( sol.synthesisDL.D * sol.synthesisDL.D' ); 
% % tmp1 = U2 * sol.synthesisDL.D;
% 
% 
% 
% 
% 
% figure; imagesc( abs(synthesisDL.D' * synthesisDL.D ))
% 
% Up = padarray( U, [0, sol.synthesisDL.Na - sol.synthesisDL.Ni], 'post' );
% synthesisDL.D = Up * V';
% figure; imagesc( abs(synthesisDL.D' * synthesisDL.D ))
% 
% 
% [U,S,V] = svd( synthesisDL.D', 0 );
% synthesisDL.D = (U * V').';
% 
% 
% 
% 
% % rescale columns ( synthesis dictionary terms ) to be unit norm:
% synthesisDL.D = synthesisDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( synthesisDL.D ) .^ 2 )), [ sol.synthesisDL.Ni, 1 ] ));
% 
% 
% [U,S,V] = svd( synthesisDL.D, 0 );
% Up = padarray( U, [0, synthesisDL.Na - synthesisDL.Ni], 'post' );
% synthesisDL.D = Up * V';
% 
% figure; imagesc( abs(synthesisDL.D' * synthesisDL.D ))
% 
% 
% 
% 
% % orthogonal/orthonormal dictionary columns:
% [ U, E ] = eig( synthesisDL.D' * synthesisDL.D ); 
% %Ds = Ds * U ./ ( 1e-6 + repmat( sqrt( diag( E )).' , [ synthesisDL.Ni, 1 ] ) ); 
% synthesisDL.D = synthesisDL.D * U; 
% 
% 
% [ U, E ] = eig( sol.synthesisDL.D' * sol.synthesisDL.D ); 
% synthesisDL.D = sol.synthesisDL.D * U; 
% figure; imagesc( abs(synthesisDL.D * synthesisDL.D' ))
% figure; imagesc( log10(1+abs(synthesisDL.D' * synthesisDL.D )))
% 
% [ U, E ] = eig( sol.synthesisDL.D' * sol.synthesisDL.D ); 
% synthesisDL.D = (sol.synthesisDL.D' * U)'; 
% figure; imagesc( abs(synthesisDL.D' * synthesisDL.D ))
% 
% return










