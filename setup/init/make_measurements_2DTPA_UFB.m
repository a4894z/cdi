function [ meas, csysz3 ] = make_measurements_2DTPA_UFB( expt )

% make_2DTPA_UFB_measurements
% make_2DRMS_measurements
% make_2DRPA_measurements
% make_2DTPA_measurements
% make_3DlaminoTPA_measurements
% make_3DlaminoTMS_measurements
% make_3DtomoTPA_measurements
% make_3DtomoTMS_measurements
% make_3DBPA_measurements

%====================================================================================================================================================
%
% define the exitwave at sample plane z2 ( assuming projection 
% approximation ) and propagate from z2 to measurement plane z3 
% using the "direct" method for long range wave field propagation
%
%============================================================
% measurement plane z3 real/reciprocal space pixel sizes, etc 
%============================================================

% fov of measurement plane z3 using sample plane z2 parameters:
csysz3.Lr = single( expt.lambda * expt.csys.z23 / expt.csys.z2.dLy ); 
csysz3.Lc = single( expt.lambda * expt.csys.z23 / expt.csys.z2.dLx ); 

% pixel sizes at the focal plane z1:
csysz3.dr = csysz3.Lr / expt.sz.r;
csysz3.dc = csysz3.Lc / expt.sz.c;

% max spatial frequencies for detector:
csysz3.qrmax = 1 / csysz3.dr;
csysz3.qcmax = 1 / csysz3.dc;

% spatial frequency pixel size:
csysz3.dqr = csysz3.qrmax / expt.sz.r;
csysz3.dqc = csysz3.qcmax / expt.sz.c;
 
%======================
% sample in measurement 
%======================

% figure; 
% plot_2Dscan_positions( expt.spos.rs, [], expt.spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-500, 500])
% ylim([-500, 500])
% daspect([1 1 1])  
% grid on

for ss = 1 : expt.spos.N
    
    fprintf( [ num2str( [ ss, expt.spos.N ], 'Creating measurement for scan position %d / %d' ), '\n' ]);
    
    % compute exitwave(s) for the SCPMs and sample at the current scan position:
    [ phiz2, ~ ] = enforce_2DTPAsposview( expt.probe.phi, expt.sample.T, expt.sample.vs.r, expt.sample.vs.c, expt.spos.rs( ss, : ), 'subpx' );
    
    [ phiz3 ] = wfprop2DTPA_UFB_z2z3( phiz2, expt );

    img = sum( abs( phiz3 ) .^ 2, 3 );
   
    % intensity --> magnitude
    sqrt_img = sqrt( img );

    % compute norms of the measurements
    meas.SI_sumD2( ss )         = sum( img( : ));          
    meas.SI_fronormD( ss )      = sqrt( meas.SI_sumD2( ss ));      % norm( sqrt_img, 'fro' );
    
    meas.D( :, :, ss )     =  sqrt_img;
    meas.Deq0( :, :, ss )  = ( meas.D( :, :, ss ) == 0 ); 
end

clearvars -except expt meas csysz3
fprintf('\n');

%================================
% introduce q = 0 shifting errors
%================================

% for ss = 1 : expt.spos.N
%     
%     meas.SI( ss ).D = fftshift( meas.SI( ss ).D );
%     meas.SI( ss ).D = zeropadarray( meas.SI( ss ).D, [ 5, 5 ] );
%     
% %     meas.SI( ss ).D = nocircshift2D( meas.SI( ss ).D, round( 1 * [ 2 * rand - 1, 2 * rand - 1 ] ));
%     meas.SI( ss ).D = nocircshift2D( meas.SI( ss ).D, round( 1 * [ 2, -1 ] ));
% %     meas.SI( ss ).D = subpixelshift2D( meas.SI( ss ).D, [1.5, 0] );
% 
%     meas.SI( ss ).D = truncatearray( meas.SI( ss ).D, expt.sz.sz );
%     meas.SI( ss ).D = fftshift( meas.SI( ss ).D );
%     
% %     figure( 666 ); 
% %     imagesc( log10(1+meas.SI( ss ).D) )
% %     pause( 0.1 )
%     
% end

%========================================
% beamstop, dead pixels, attenuation, etc
%========================================

% for ss = 1 : expt.spos.N
% 
%     tmp0 = fftshift( meas.SI( ss ).D );
%     tmp0( end - 40 : end, : ) = 0;  
% %     tmp0( 1 : 40, : ) = 0;  
%     meas.SI( ss ).D = fftshift( tmp0 );
%     
%     
%     
%     
%     figure( 666 ); 
%     imagesc( fftshift( log10(1+meas.SI( ss ).D) ));
%     daspect([1 1 1]);
%     colormap jet
%     pause( 0.05 )
%     
%     
% 
% end

%====================================
% introduce noise to the measurements
%====================================

% [ snr_reconst, phot_px ] = compute_snr_vs_q( fftshift( meas.D( :, :, 353 )) );
% figure( 666 ); 
% plot( phot_px )

meas.noisy = 1;

if logical( meas.noisy ) == true

    meas.Nexposures = 3;
    
    %=======================
    % sample in measurements 
    %=======================

    meas_noisy = meas.D .^ 2;
    csum       = 0;
    sz         = size( meas.D );
    max_BG     = 1;                    % define max background value allowed ( for unif random distribution ):

    %========
    
    meas_noisy = gpuArray( meas_noisy );
    csum       = gpuArray( csum );
    sz         = gpuArray( sz );
    max_BG     = gpuArray( max_BG );

    %========
    
    for ii = 1 : meas.Nexposures 
        
        fprintf( [ num2str( [ ii, meas.Nexposures ], 'Contaminating Measurements with Noise, Exposure = %d / %d' ), '\n' ]);
        
        csum = csum + random( 'poisson', meas_noisy ) + random( 'unif', 0, max_BG, sz );

    end

    meas_noisy = gather( csum / meas.Nexposures );

    %===================
    % Background removal
    %===================
    
    bg_sub = 0.95 * max_BG;   
    
%     meas_noisy = meas_noisy - bg_sub;
%     meas_noisy( meas_noisy < 0 ) = 0;

%     meas_noisy( meas_noisy < 0.95e-4 * max( meas_noisy(:) ) ) = 0;
    meas_noisy( meas_noisy < bg_sub ) = 0;
    
    for ss = 1 : 3 : size( meas.D, 3 )
        
        figure( 666 ); 
        subplot(121); imagesc(log10( 1 + 10^0 * fftshift(abs(meas.D(:,:,ss) .^ 2)))); daspect([1 1 1])
        subplot(122); imagesc(log10( 1 + 10^0 * fftshift(abs(meas_noisy(:,:,ss))))); daspect([1 1 1])
        colormap turbo
        
    end
%     norm(  meas.D(:)  ) / norm( sqrt( meas_noisy(:) ) )

    %========
        
    meas.D    = gather( sqrt( meas_noisy ));
    meas.Deq0 = ( meas.D == 0 );
 
end

% [ snr_reconst, phot_px ] = compute_snr_vs_q( fftshift( meas.D( :, :, 353 )) );
% figure( 666 ); 
% hold on
% plot( phot_px )
% hold off
% 
% 5;

