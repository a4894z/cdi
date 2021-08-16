%

%{

clear; close all; 
restoredefaultpath
codepath = '~/Documents/MATLAB/Code/cdi/';   
cd( codepath );
addpath(genpath(pwd));
clear


clear; close all; testing_reflectiongeodata_s8

%}

%==================================================================================================
%-------------------------------- SET UP PATHS, LOAD DATA, ETC ------------------------------------
%==================================================================================================

fprintf('\n========================================================================================================'); 
fprintf('\nLOADING Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed...'); 

load( '~/Documents/MATLAB/Code/cdi/remap_Surface_A4_star_phi001_0.0.mat' );

fprintf('done loading data!\n'); 
fprintf('========================================================================================================\n\n'); 

% expt.paths.code = codelocation;
% expt.paths.rrundata = rrundata; 
% expt.paths.srundata = rrundata; 

% % check if results directory exists; if not make it
% sol.paths.imresults = [ codelocation, 'imresults' ];
% if exist( sol.paths.imresults, 'dir' ) ~= 7, mkdir( sol.paths.imresults ), end

% clearvars -except expt sol

%==================================================================================================

sol.plot.xaxis = expt.csys.z2.dx * (1 : sol.sz.c) * 1e6;
sol.plot.xaxis = sol.plot.xaxis - mean( sol.plot.xaxis );
sol.plot.yaxis = expt.csys.z2.dy * (1 : sol.sz.r) * 1e6; 
sol.plot.yaxis = sol.plot.yaxis - mean( sol.plot.yaxis );
 
sol.sparse.phiedges = edge_detect_xy_forward_diff( [], sol.sz );
sol.probe.scpm.N = 1;

tmp0 = make_gaussian( sol.sz.sz, [ 0.5 * sol.sz.r + 1, 0.5 * sol.sz.c + 1 ], [ 0.5 * sol.sz.r, 0.5 * sol.sz.c ]);
sol.measLPF = 0 + 1 * fftshift( tmp0 );


% view slices for when computing the exit wave views at a particular scan position
% sol.sample.vs.r = round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r )));
sol.sample.vs.c = round( ( 0.5 * sol.sz.c -  25 + 1 ) :  ( 0.5 * sol.sz.c + 25 ));

%==================================================================================================

sol.ss = 47;

if ~isfield( sol, 'phi' )
    
    % sol.phi = fftshift( ifft2( meas.D )) * expt.sz.sqrt_rc;
    % sol.phi = make_gaussian( sol.sz.sz, [ 0.5 * sol.sz.r + 1, 0.5 * sol.sz.c + 1 ], [ 0.5 * sol.sz.r, 0.5 * sol.sz.c ]);
    sol.phi = rand( sol.sz.sz, 'single' ) + 1i * rand( sol.sz.sz, 'single' );
    
end

% CENTER THE EXIT WAVE USING MAX VALUE:
[ cr, cc ] = find( sol.phi == max( max( sol.phi )));
sol.phi = circshift( sol.phi, [128, 256] - [ cr, cc ] );

% lagrangian multipliers
lam.x = complex( zeros( expt.sz.sz, 'single' ));
lam.y = complex( zeros( expt.sz.sz, 'single' ));

[ Xm ] = forward_edgedetect( sol.phi, sol.sparse.phiedges );

%==================================================================================================

% TESTING: SPARSE EXIT WAVE, SINGLE VIEW:

for ii = 1 : 1000
  
    fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
    if ( mod( sol.it.exwv, 25 ) == 0 ), fprintf( '\n' ); end
        
    %==============================================================================================

    Xs.x = Xm.x - lam.x;
    Xs.y = Xm.y - lam.y;
    
    %==============================================================================================
    % ISO ABS SPARSITY
    
    %%{
    
    abs_Xs.x = abs( Xs.x );
    abs_Xs.y = abs( Xs.y );
    
    abs_x2y2 = sqrt( abs_Xs.x .^ 2 + abs_Xs.y .^ 2 );
    
    if mod( ii, 1e9 ) == 0

        sol.sparse.philvl = round( sol.sz.rc * 0.92 );
        lamx2y2 = find_thresh_from_sparsitylevel( abs_x2y2, sol.sparse.philvl );
        
        [ Xs.x, Wr ] = soft_shrinkage( Xs.x, abs_x2y2, lamx2y2 );
        [ Xs.y, Wc ] = soft_shrinkage( Xs.y, abs_x2y2, lamx2y2 );


    else

        sol.sparse.philvl = round( sol.sz.rc * 0.07 );
        lamx2y2 = find_thresh_from_sparsitylevel( abs_x2y2, sol.sparse.philvl );

        [ Xs.x, Wc ] = hard_shrinkage( Xs.x, abs_x2y2, lamx2y2 );
        [ Xs.y, Wr ] = hard_shrinkage( Xs.y, abs_x2y2, lamx2y2 );


    end
    
    %}
    
    %==============================================================================================
    % ANISO ABS SPARSITY
    
    %{
    
    abs_Xs.x = abs( Xs.x );
    abs_Xs.y = abs( Xs.y );

    if mod( ii, 5 ) == 0

        sol.sparse.philvl = round( sol.sz.rc * 0.3 );
        lamc = find_thresh_from_sparsitylevel( abs_Xs.x, sol.sparse.philvl );
        lamr = find_thresh_from_sparsitylevel( abs_Xs.y, sol.sparse.philvl );

        [ Xs.x, Wr ] = soft_shrinkage( Xs.x, abs_Xs.x, lamc );
        [ Xs.y, Wc ] = soft_shrinkage( Xs.y, abs_Xs.y, lamr );


    else

        sol.sparse.philvl = round( sol.sz.rc * 0.15 );
        lamc = find_thresh_from_sparsitylevel( abs_Xs.x, sol.sparse.philvl );
        lamr = find_thresh_from_sparsitylevel( abs_Xs.y, sol.sparse.philvl );

        [ Xs.x, Wc ] = hard_shrinkage( Xs.x, abs_Xs.x, lamc );
        [ Xs.y, Wr ] = hard_shrinkage( Xs.y, abs_Xs.y, lamr );


    end
    
    %}
    
    %==============================================================================================
    
    %{
    
    re_Xs.x = real( Xs.x );
    re_Xs.y = real( Xs.y );
    im_Xs.x = imag( Xs.x );
    im_Xs.y = imag( Xs.y );

    abs_re_Xs.x = abs( re_Xs.x );
    abs_re_Xs.y = abs( re_Xs.y );
    abs_im_Xs.x = abs( im_Xs.x );
    abs_im_Xs.y = abs( im_Xs.y );

    if mod( ii, 1 ) == 0

        sol.sparse.philvl = round( sol.sz.rc * 0.43 );
        lam_rex = find_thresh_from_sparsitylevel( abs_re_Xs.x, sol.sparse.philvl );
        lam_rey = find_thresh_from_sparsitylevel( abs_re_Xs.y, sol.sparse.philvl );
        lam_imx = find_thresh_from_sparsitylevel( abs_im_Xs.x, sol.sparse.philvl );
        lam_imy = find_thresh_from_sparsitylevel( abs_im_Xs.y, sol.sparse.philvl );

        [ re_Xs.x, Wrex ] = soft_shrinkage( re_Xs.x, abs_re_Xs.x, lam_rex );
        [ re_Xs.y, Wrey ] = soft_shrinkage( re_Xs.y, abs_re_Xs.y, lam_rey );
        [ im_Xs.x, Wimx ] = soft_shrinkage( im_Xs.x, abs_im_Xs.x, lam_imx );
        [ im_Xs.y, Wimy ] = soft_shrinkage( im_Xs.y, abs_im_Xs.y, lam_imy );

    else

        sol.sparse.philvl = round( sol.sz.rc * 0.05 );
        lam_rex = find_thresh_from_sparsitylevel( abs_re_Xs.x, sol.sparse.philvl );
        lam_rey = find_thresh_from_sparsitylevel( abs_re_Xs.y, sol.sparse.philvl );
        lam_imx = find_thresh_from_sparsitylevel( abs_im_Xs.x, sol.sparse.philvl );
        lam_imy = find_thresh_from_sparsitylevel( abs_im_Xs.y, sol.sparse.philvl );

        [ re_Xs.x, Wrex ] = hard_shrinkage( re_Xs.x, abs_re_Xs.x, lam_rex );
        [ re_Xs.y, Wrey ] = hard_shrinkage( re_Xs.y, abs_re_Xs.y, lam_rey );
        [ im_Xs.x, Wimx ] = hard_shrinkage( im_Xs.x, abs_im_Xs.x, lam_imx );
        [ im_Xs.y, Wimy ] = hard_shrinkage( im_Xs.y, abs_im_Xs.y, lam_imy );


    end

    Xs.x =  re_Xs.x + 1i * im_Xs.x;
    Xs.y =  re_Xs.y + 1i * im_Xs.y;
    
    %}
    
    %==============================================================================================

    Xm.x = Xs.x + lam.x;
    Xm.y = Xs.y + lam.y;

    sol.phi = inverse_edgedetect( Xm, sol.sparse.phiedges, sol.sz );

    sol.phi = enforce_2DTPAmeas( sol.phi, expt.meas.SI( sol.ss ), sol.measLPF, sol );

    %===============

    [ Xm ] = forward_edgedetect( sol.phi, sol.sparse.phiedges );

    %===============

    lam.x = lam.x + 1.0 * ( Xs.x - Xm.x );
    lam.y = lam.y + 1.0 * ( Xs.y - Xm.y );

    %===============

    if mod( sol.it.exwv, 25 ) == 0 

        figure( 666 ); 
        set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        subaxis(2,1,1,'SpacingVert',0,'MR',0.02, 'ML',0.04,'MT',0.03, 'MB',0.13); 
        imagesc( sol.plot.xaxis, sol.plot.yaxis, log10(1+abs(sol.phi))); 
        xlabel('um'); ylabel('um'); grid on; colormap(expt.cm.blj)
        subaxis(2,1,2,'SpacingVert',0,'MR',0.02, 'ML',0.04,'MT',0.13, 'MB',0.05);
        imagescHSV(log10(1+abs(sol.phi)) .* exp(1i*angle(sol.phi)), sol.plot ); 
        xlabel('um'); ylabel('um'); grid on
        title(num2str(sol.it.exwv,'%d'))
        
        export_fig( num2str( [ sol.it.exwv, sol.ss ], 'exitwave_%d_spos-%d.jpg' ), '-r90.0' )
        close all;
        
    end

    %==================================
    % EXITWAVE ITERATION COUNTER
    %==================================
    
    sol.it.exwv = sol.it.exwv + 1;
    
end

%==================================================================================================
%---------------------------------- SHRINKWRAP + HIO FOR EXITWAVE ---------------------------------
%==================================================================================================

for ii = 1 : 0
    
    %ii
    
    phiM = enforce_2DTPAmeas( sol.phi, expt.meas.SI( sol.ss ), sol.measLPF, sol );

    abs_phi = abs( sol.phi );
    blur_abs_phi = lpf_gauss( abs_phi, 0.02 * sol.sz.sz );
    blur_abs_phi = abs( blur_abs_phi );

    thresh = find_thresh_from_sparsitylevel( blur_abs_phi, round( sol.sz.rc * 0.025 ) );
    W = ( blur_abs_phi > thresh );
    
    sol.phi = W .* phiM;
    
%     if mod( ii, 50e9 ) ~= 0 
%         
%         % HIO
%         sol.phi = W .* phiM + not( W ) .* ( sol.phi - 0.8 * phiM );
%         
%     else
%         
%         % ER
%         sol.phi = W .* phiM;
%         
%     end

    % sol.phi = pis pim sol.phi + not(pis) (sol.phi - beta * pim sol.phi )

    
    figure(666); 
    subplot(121); imagesc(log10(1+abs(sol.phi))); daspect([1 1 1])
    subplot(122); imagescHSV(log10(1+abs(sol.phi)) .* exp(1i*angle(sol.phi))); daspect([1 1 1])
    title(num2str(ii,'%d'))
    pause( 0.05)
    

end

%==================================================================================================

clearvars -except expt sol

% sol.phi = []; sol = rmfield( sol, 'phi' );
% sol.phiOLD = []; sol = rmfield( sol, 'phiOLD' );
% sol.unmeas = []; sol = rmfield( sol, 'unmeas' );

% '~/Documents/MATLAB/Code/cdi/remap_Surface_A4_star_phi001_0.0.mat'


% save data
fprintf('\n========================================================================================================'); 
fprintf('\nSAVING Simulated Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed...'); 
save( expt.paths.sdata, '*' );
% save( [ expt.paths.sdata, sprintf('Star_fine_mid_out_%s.mat', datestr(now,'ddmmmyyyy_HHMMSS')) ], '*' );
% probe = sol.probe; 
% sample = sol.sample; 
% spos = sol.spos; 
% save( [ expt.paths.sdata, num2str( sol.it.exwv,' Star_fine_mid_out_it%d '), sprintf('_%s.mat', datestr(now,'ddmmmyyyy_HHMMSS')) ], 'probe', 'sample', 'spos' );
fprintf('done saving data!\n'); 
fprintf('========================================================================================================\n\n'); 

