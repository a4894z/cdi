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

%==================================================================================================
%
% DESCRIPTION:
%
% define the exitwave at sample plane z2 ( assuming projection 
% approximation ) and propagate from z2 to measurement plane z3 
% using the "direct" method for long range wave field propagation
%
%==================================================================================================
%-------------- measurement plane z3 real/reciprocal space pixel sizes, etc -----------------------
%==================================================================================================

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

%==================================================================================================
%-------------------------------------- sample out measurement ------------------------------------
%==================================================================================================
 
% wfp.lambda = expt.lambda; 
% wfp.zif = expt.csys.z23; 
% wfp.zi = expt.csys.z2;
% wfp.zf = expt.csys.z3;
% wfp.intcurv = true;
% wfp.extcurv = false;
% wfp.dir = 'forward';
% tmp1 = fresnelpropdirect( expt.probe.phi, wfp );

% % TODO: USE FUNCTION HANDLE HERE!!!!!!!!!!!!!!!!!!!!!!
% % propagate probe modes from sample plane z2 to measurement plane z3:
% [ probez3 ] = wfprop2DTPA_UFB_z2z3( expt.probe.phi, expt );
% 
% cumsum = zeros( expt.sz.r, expt.sz.c, 'single');
% 
% % probe mode wavefield intensities at z3 sum incoherently: 
% for pp = 1 : expt.probe.scpm.N
% 
%     cumsum = cumsum + abs( probez3( :, :, pp )) .^ 2;
%     
% %     figure; 
% %     subplot(131); imagesc_diffraction( cumsum ); colormap( expt.cm.blj )
% %     subplot(132); imagesc_diffraction( probez3( :, :, pp ) ); colormap( expt.cm.blj )
% %     subplot(133); imagescHSV(  expt.probe.phi( :, :, pp )); daspect([1 1 1])
%     
% end
%                         
% meas.SO = sqrt( cumsum );



    


    
%     img = fft2( fftshift( expt.sample.T )) / expt.sample.sz.sqrt_rc;
%     norm_img = norm( img, 'fro' ); 
%     
%     [ xx, yy ] = meshgrid( 1 : expt.sample.sz.c, 1 : expt.sample.sz.r );
%     qintensityscaling = 1 * yy.^4 + 0 * xx.^3;
%     
%     [ qintensityscaling ] = modulus_limits_scale( qintensityscaling, [ 0.05, 8.0 ] );
%     
%     figure; 
%     imagesc( fftshift( log10(1+ 10^2 * abs( img )))); 
%     daspect([1 1 1]); 
%     colormap(jet);
%     colorbar;
%     
%     img = img .* fftshift( qintensityscaling );
%     img = img * norm_img / norm( img, 'fro' );
%     
%     figure; 
%     imagesc( fftshift( log10(1+ 10^2 * abs( img )))); 
%     daspect([1 1 1]); 
%     colormap(jet);
%     colorbar;
%    
%     sample_TF_qscaling = fftshift( ifft2( img )) * expt.sample.sz.sqrt_rc;
%     
%     
% 
%    
%     figure; 
%     imagesc( abs( expt.sample.T )); daspect([1 1 1]); 
%     colormap(jet)
%     colorbar;
%     
%     figure; 
%     imagescHSV( ( sample_TF_qscaling )); daspect([1 1 1]); 
%     colormap(jet)
%     colorbar;
%     
%     figure; imagesc( qintensityscaling ); 
%     daspect([1 1 1]); 
%     colormap(jet)
%     colorbar;
% 
%     figure; 
%     imagesc( fftshift( abs(ifft2( qintensityscaling ) * expt.sample.sz.sqrt_rc ))); 
%     daspect([1 1 1]); 
%     colormap(jet)
%     colorbar;
    
    
    
%==================================================================================================
%-------------------------------------- sample in measurement -------------------------------------
%==================================================================================================

for ss = 1 : expt.spos.N
    
    fprintf( [ num2str( [ ss, expt.spos.N ], 'Creating measurement for scan position %d / %d' ), '\n' ]);
    
    % compute exitwave(s) for the SCPMs and sample at the current scan position:
    [ phiz2, ~ ] = enforce_2DTPAsposview( expt.probe.phi, expt.sample.T, expt.sample.vs.r, expt.sample.vs.c, expt.spos.rs( ss, : ), 'subpx' );
    
    % TODO: USE FUNCTION HANDLE HERE!!!!!!!!!!!!!!!!!!!!!!
    % numerically propagate these exitwaves from the sample plane z2 to the measurement plane z3:
    [ phiz3 ] = wfprop2DTPA_UFB_z2z3( phiz2, expt );
    % wfprop_z2z3UFB_2DTPA

    
    
    

    
    
    
    
    
    
    
    img = sum( abs( phiz3 ) .^ 2, 3 );
    
%     figure( 666 );
%     imagesc( log10(1 + fftshift( img )));
%     daspect([1 1 1])
%     colormap gray
    
    %========================
%     
%     norm_img = norm( img, 'fro' ); 
%     
% %     [ xx, yy ] = meshgrid( ( 1 - 0.5 * expt.sz.c ) : 0.5 * expt.sz.c, ( 1 - 0.5 * expt.sz.r ) : 0.5 * expt.sz.r );
%     [ xx, yy ] = meshgrid( 1 : expt.sz.c, 1 : expt.sz.r );
%     
%     qintensityscaling = yy.^4;
%     
%     [ qintensityscaling ] = modulus_limits_scale( qintensityscaling, [ 0.05, 3.0 ] );
%     
% %     qintensityscaling = qintensityscaling / max( abs( qintensityscaling( : )));
% %     qintensityscaling = qintensityscaling + 0.3;
% %     qintensityscaling = qintensityscaling / max( abs( qintensityscaling( : )));
%     
%     img = img .* fftshift( qintensityscaling );
%     
%     img = img * norm_img / norm( img, 'fro' );
    
    %========================
    
%     % GET RID OF FOR LOOP
%     cumsum = zeros( expt.sz.r, expt.sz.c, 'single');
%         
%     % probe mode wavefield intensities at z3 sum incoherently: 
%     for pp = 1 : expt.probe.scpm.N
% 
%         cumsum = cumsum + abs( phiz3( :, :, pp )) .^ 2;
%         
% %         figure; 
% %         imagesc( fftshift( log10( 1 + 10^2 * cumsum )), [ 0, 7 ])
% %         daspect([1 1 1])
% %         colormap( expt.cm.blj )
%         
%     end
    
    



    
%     img = nocircshift2D( img, round( [ 2 * ( 2 * rand - 1), 2 * ( 2 * rand - 1) ]));
    
    
    
    
    
    
    
    % intensity --> magnitude
    sqrt_img = sqrt( img );

    % compute norms of the measurements
    meas.SI_sumD2( ss )         = sum( img( : ));          
    meas.SI_fronormD( ss )      = sqrt( meas.SI_sumD2( ss ));      % norm( sqrt_img, 'fro' );

%     % assign the measurements to relevant arrays, record missing diffraction data
%     meas.SI( ss ).D            = sqrt_img;
%     meas.SI( ss ).Deq0         = ( meas.SI( ss ).D == 0 );
%     meas.SI( ss ).Dneq0        = not( meas.SI( ss ).Deq0 );

    
    
    
    
    
    
        meas.D( :, :, ss )     =  sqrt_img;
        meas.Deq0( :, :, ss )  = ( meas.D( :, :, ss ) == 0 );
%         meas.Dneq0( :, :, ss ) = not( meas.Deq0( :, :, ss ) );
        
        
        
        
        
end

clearvars -except expt meas csysz3
fprintf('\n');

%==================================================================================================
%---------------------------- introduce q = 0 shifting errors -------------------------------------
%==================================================================================================

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

%==================================================================================================
%----------------------- beamstop, dead pixels, attenuation, etc ----------------------------------
%==================================================================================================

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

%==================================================================================================
%------------------------- introduce noise to the measurements ------------------------------------
%==================================================================================================

% measurments just defined are noise free, now introduce poisson noise + detector background noise, 
% i.e. a noise model like: In = Poiss[ I ] + Unif[ ]

meas.noisy = 1;

if logical( meas.noisy ) == true

    % total number of exposures for sample in measurement:
    meas.NexposuresSI = 1;

    % total number of exposures for sample out measurement:
    meas.NexposuresSO = 1;
 
    %=======================
    % sample in measurements 
    %=======================

    meas_noisy = meas.D .^ 2;
    csum       = 0;
    sz         = size( meas.D );
    max_BG     = 10;                    % define max background value allowed ( for unif random distribution ):
    Nchunk     = 200;

    %============================================
    
    meas_noisy = gpuArray( meas_noisy );
    csum       = gpuArray( csum );
    sz         = gpuArray( sz );
    max_BG     = gpuArray( max_BG );
    Nchunk     = gpuArray( Nchunk );
    
    %============================================
    
    for ii = 1 : meas.NexposuresSI 
        
        fprintf( [ num2str( [ ii, meas.NexposuresSI ], 'Contaminating Measurements with Noise, Exposure = %d / %d' ), '\n' ]);
        
        csum = csum + random( 'poisson', meas_noisy ) + random( 'unif', 0, max_BG, sz );

    end

    meas_noisy = gather( csum / meas.NexposuresSI );
    
    %============================================
    
%     Nchonk = gpuArray( floor( sz( 3 ) /  Nchunk ));
% 
%     for jj = 1 : Nchonk
% 
%         chonk = ( ( jj - 1 ) * Nchunk + 1 ) : ( jj * Nchunk );
% 
%         for ii = 1 : meas.NexposuresSI 
% 
%             csum = csum + random( 'poisson', meas_noisy( :, :, chonk )) + random( 'unif', 0, max_BG, [ sz( 1 ), sz( 2 ), Nchunk ] );
% 
%         end
% 
%         meas_noisy( :, :, chonk ) = csum;
% 
%         csum = gpuArray( 0 );
% 
%     end
%     
%     %==============
%     
%     chonk  = ( 1 + chonk( end )) : sz( 3 );
%     Nchunk = length( chonk );
% 
%     for ii = 1 : meas.NexposuresSI 
% 
%         csum = csum + random( 'poisson', meas_noisy( :, :, chonk )) + random( 'unif', 0, max_BG, [ sz( 1 ), sz( 2 ), Nchunk ] );
% 
%     end
% 
%     meas_noisy( :, :, chonk ) = csum;
%     
%     %==============
% 
%     meas_noisy = gather( meas_noisy / meas.NexposuresSI );
% 
%     clear( 'chonk', 'csum', 'sz', 'max_BG', 'Nchunk', 'Nchonk', 'ii', 'jj' )

    %===================
    % Background removal
    %===================
    
%     bg_sub = 1 * 4.0;                       % use 4.0 for x 10^-1 wrt orig
    bg_sub = 0 * 8.0;                       % use 8.0 for x 10^+0 wrt orig
%     bg_sub = 0 * 8.0;                       % use 8.0 for x 10^+1 wrt orig

    meas_noisy = meas_noisy - bg_sub;
    
    %==============
    
    meas_noisy( meas_noisy < 0 ) = 0;

    %==============
    
    figure; imagesc(log10( 1 + 10^0 * fftshift(abs(meas.D(:,:,155) .^ 2))))
    
    figure; imagesc(log10( 1 + 10^0 * fftshift(abs(meas_noisy(:,:,155)))))
    
    5;
    
    meas.D    = sqrt( meas_noisy );
    meas.Deq0 = meas_noisy == 0;

%     figure; imagesc(log10( 1 + 10^2 * fftshift(abs(meas.D(:,:,195) .^ 2 ))))

    
    
    %=======================
    % Sample OUT measurement 
    %=======================

%     % magnitude --> intensity
%     meas_SO = single( meas.SO .^ 2 );
% 
%     tmpSO = zeros( expt.sz.r, expt.sz.c, 'single' );
% 
%     for ii = 1 : meas.NexposuresSO
% 
%         % TODO: add in an effect where the wavefield before the focusing optic has some variation shot to shot
%         % e.g. a blurred random array with variations of +/- 5% or so?
%         % or a gaussian of some fwhm that can "jitter" spatially wrt optic position
% 
%         % poisson noise in intensity measurements:
%         tmp2 = random( 'poisson', meas_SO );
% 
%         % constant noisy background ( uniform random distributed ):
%         tmp1 = random( 'unif', 0, max_BG, expt.sz.r, expt.sz.c );
% 
%         % accumulate single exposure intensity measurements:
%         tmpSO = tmpSO + tmp2 + tmp1;
% 
%     end
% 
%     meas.SO = sqrt( tmpSO / meas.NexposuresSO );



















%     % magnitude --> intensity
%     meas.D = meas.D .^ 2;
% 
%     tmpSI = zeros( expt.sz.r, expt.sz.c, 'single' );
%     
%     for ss = 1 : expt.spos.N
%         
%         fprintf( [ num2str( [ ss, expt.spos.N, meas.NexposuresSI ], 'Introducing measurement noise for scan position %d / %d, No. of Exposures = %d' ), '\n' ]);
%         
%         for ii = 1 : meas.NexposuresSI  
% 
%             % TODO: add in an effect where the wavefield before the focusing optic has some variation shot to shot
%             % e.g. a blurred random array with variations of +/- 5% or so?
%             % or a gaussian of some fwhm that can "jitter" spatially wrt optic position
% 
%             % poisson noise in intensity measurements:
%             tmp1 = random( 'poisson', meas.D( :, :, ss ) );
% 
%             % constant noisy background ( uniform random distributed ):
%             tmp2 = random( 'unif', 0, max_BG, expt.sz.r, expt.sz.c );
% 
%             % accumulate single exposure intensity measurements:
%             tmpSI = tmpSI + tmp1 + tmp2;
% 
%         end
% 
%         meas.D( :, :, ss ) = sqrt( tmpSI / meas.NexposuresSI );
%         tmpSI = zeros( expt.sz.r, expt.sz.c, 'single' );
%     
%         figure(555); 
%         imagesc_diffraction( meas.D( :, :, ss )); daspect([1 1 1]); colormap( expt.cm.blj )
%         drawnow
%         
% %         close all
%         
%     end
%     
%     fprintf('\n');
    
end

%===================
% Background removal
%===================

% meas.rmbackgrnd = 0;
% 
% if logical( meas.rmbackgrnd ) == true
%     
%     % constant background subtraction on intensity:
%     bg_sub = 8.0;
%     
%     for ss = 1 : expt.spos.N
%         
%         tmp0 = abs( meas.D( :, :, ss ) ) .^ 2 - bg_sub;
%         tmp0( tmp0 < 0 ) = 0;
%         tmp0 = sqrt( tmp0 + 0 * ( tmp0 ~= 0 ) * bg_sub );
%      
%         figure( 666 ); 
%         set( gcf, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
%         subplot(121);
%         imagesc( fftshift( log10( 1 + meas.D( :, :, ss ) ))); daspect([1 1 1]); colormap( expt.cm.blj )
%         daspect([1 1 1])
%         subplot(122);
%         imagesc( fftshift( log10( 1 + tmp0 ))); daspect([1 1 1]); colormap( expt.cm.blj )
%         daspect([1 1 1])
%         pause( 0.06 )
%         
%         meas.D( :, :, ss ) = tmp0;
%         
%     end
%        
% %         figure; 
% %         imagesc_diffraction( meas.SO ); daspect([1 1 1]); colormap( expt.cm.blj )
%         
%     bg_sub = 10.0;
%     meas.SO = abs( meas.SO ) .^ 2 - bg_sub;
%     meas.SO( meas.SO < 0 ) = 0;
%     meas.SO = sqrt( meas.SO + 0 * (meas.SO ~= 0) * bg_sub );
%     
% %         figure(222); 
% %         imagesc_diffraction( meas.SO ); daspect([1 1 1]); colormap( expt.cm.blj )
%       
% end

%================================================

% % missing data in measurement(s):
% 
% meas.SOeq0 = ( meas.SO == 0 );
% meas.SOneq0 = not( meas.SOeq0 );
% 
% for ss = 1 : expt.spos.N
%     meas.SI( ss ).Deq0 = ( meas.D( :, :, ss ) == 0 );
%     meas.SI( ss ).Dneq0 = not( meas.SI( ss ).Deq0 );
% end

%================================================

% % clean up 
% clearvars -except expt meas
