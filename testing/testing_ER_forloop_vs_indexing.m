%

%{

clear; close all; 

restoredefaultpath
codelocation =  '~/Documents/Science/Matlab/Code/cdi/';
cd( codelocation );
addpath( genpath( pwd ));   
clearvars -except expt sol



clear; close all; testing_ER_forloop_vs_indexing

%}

%==================================================================================================

clear; close all;
load /home/ash/Documents/Science/Matlab/Code/cdi/sim_ptycho2DTPA.mat;

[ sol, expt ] = set_defaults( sol, expt );  

sol.probe  = expt.probe;
% sol.sample = expt.sample;
% sol.spos   = expt.spos;
% sol.sz     = expt.sz;

clearvars -except expt sol

sol.spos.updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));
sol.spos.shifttype   = 'px';

%================================================



    


    
    
    
    
    
    
    
    
    
tic 
    
 
[ sol.phi] = ERupdate_2DTPAexwv_PTFspos_meas( sol, expt );

  
[ sample1 ] = DMupdate_2Dsample_loop( sol );

[ P ] = DMupdate_probemodes( sample1, sol.phi, sol );
toc

% figure; imagesc( abs(sample1), [0, 1] )
% figure; imagesc( abs(P(:,:,3)), [0, 1] )








tic

sol.phi = ER_CPU_arrays_hadamard( sol, expt );


[ sample2 ] = DMupdate_2Dsample_loop( sol );



[ P ] = DMupdate_probemodes( sample2, sol.phi, sol );
toc





% figure; imagesc( abs(sample2), [0, 1] )
figure; imagesc( abs(P(:,:,3)), [0, 44] )
figure; imagesc( abs(sol.probe.P(:,:,3)), [0, 44] )

5;

    %============================================
    
%     probe = gpuArray( sol.probe.P );
%     TF    = gpuArray( sol.sample.TF );
%     
%     sz = gpuArray( sol.sz.sz );
%     Nspos = gpuArray( sol.spos.N );
%     Nscpm = gpuArray( sol.probe.scpm.N );
%     
%     spos_rs = gpuArray( sol.spos.rs );
%     
%     meas_D = gpuArray( expt.meas.D );
%    
%     vs.r = gpuArray( sol.sample.vs.r );
%     vs.c = gpuArray( sol.sample.vs.c );
%     
%     
%   
% tic
% phiVg2 = ER_GPU_forloop( probe, TF, vs, spos_rs, sz, Nspos, Nscpm, meas_D );
% 
% sol.phi = gather( phiVg2 );
% 
% [ sample3 ] = DMupdate_2Dsample_loop( sol );
% 
% toc  
% 
% 
% figure; imagesc( abs( sample3 ), [0, 1] )




    %============================================
    
    sz    = gpuArray( sol.sz.sz );
    rc    = gpuArray( sol.sz.rc );
    samsz = gpuArray( sol.sample.sz.sz );
%     samrc = samsz( 1 ) * samsz( 2 );
    samrc = gpuArray( sol.sample.sz.rc );
    Nspos = gpuArray( sol.spos.N );
    Nscpm = gpuArray( sol.probe.scpm.N );
    
    ind        = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
    ind_offset = uint32( 0 : samrc : ( samrc * ( Nspos - 1 )));
%     ind_offset = 0 : samrc : ( samrc * ( spos_N - 1 ));

    probe = gpuArray( sol.probe.P );
    TFv    = gpuArray( sol.sample.TF( : ));
    TF    = gpuArray( sol.sample.TF );
    ind   = gpuArray( ind );
    
    meas_D = permute( repmat( gpuArray( expt.meas.D ), 1, 1, 1, Nscpm ), [ 1, 2, 4, 3 ] );
    
    spos_conjP_exwv = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
    spos_abs_P_abs2 = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
    

    
    
    
    
    
    


tic
phiVg = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, Nscpm, meas_D );

[ sample4b ] = DMupdate_2Dsample_repmat( spos_conjP_exwv, ...
                                         spos_abs_P_abs2, ...
                                         phiVg, ...
                                         TFv, ...
                                         probe, ...
                                         ind + ind_offset, ...
                                         rc, ...
                                         samsz, ...
                                         Nspos );



% tmp0   = sample4b( : );
tmp0   = sample4b;
TFview = reshape( tmp0( ind ), [ sz, Nspos ]);
tmp0   = permute( repmat( conj( TFview ), [ 1, 1, 1, Nscpm ] ), [ 1, 2, 4, 3 ] );
probe  = sum( tmp0 .* phiVg, 4 ) ./ sum( 1e-7 + abs( tmp0 ) .^ 2, 4 );

toc
    
    


sample4b = reshape( sample4b, sol.sample.sz.sz );


clear( 'TFv', 'ind', 'sz', 'Nspos', 'Nscpm', 'meas_D' )

sol.phi = gather( phiVg );
[ sample4 ] = DMupdate_2Dsample_loop( sol );
[ sample5 ] = DMupdate_2Dsample_loop_v2( sol );
[ sample6 ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol, expt );


figure; imagesc( abs( sample6 ), [0, 1] )
figure; imagesc( abs( sample4 ), [0, 1] )
figure; imagesc( abs( sample4b ) )
figure; imagesc( abs( sample5 ))



figure; 
subplot(221);
imagesc( abs(probe(:,:,3)), [0, 1] )
subplot(222);
imagescHSV( (probe(:,:,3)) )
subplot(223); 
imagesc( abs( sol.probe.P( :, :, 3 )), [0, 1] )
subplot(224);
imagescHSV( ( sol.probe.P( :, :, 3 )) )





% figure; imagescHSV( ( sample4 ) )
% figure; imagescHSV( ( sample4b ) )
% 
% figure; imagescHSV( (probe(:,:,3)) )
% figure; imagescHSV( (sol.probe.P(:,:,3)) )





% tic
% phiVg = ER_GPU_arrays_hadamard( p );
% toc




% ptycho_sampleupdate_DM


% phiO = ER_GPU_forloops( sol, expt );



return

%==================================================================================================
%------------------------------------------------ END ---------------------------------------------
%------------------------------- Begin utility functions for run file -----------------------------
%==================================================================================================


function [ phi ] = ER_GPU_forloop( probe, TF, vs, spos_rs, sz, Nspos, Nscpm, meas_D )

    phi = gpuArray( zeros( sz( 1 ), sz( 2 ), Nscpm, Nspos, 'single' ));      % probe .* sample view 
%     phi = zeros( sz( 1 ), sz( 2 ), Nscpm, Nspos, 'single' );      % probe .* sample view 

    sqrt_rc   = sqrt( sz( 1 ) * sz( 2 ));

    for ss = 1 : Nspos   % order * DOESN'T * matter here !!!




        %         % form 2D exitwave(s) under projection approximation in transmission geometry:
        %     [ phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( probe, TF, vs, spos_rs( ss, : ), 'px' );


        TFview = getview_2DsampleTF( TF, vs, spos_rs( ss, : ), 'px' );

        phi( :, :, :, ss ) = probe .* TFview;


        %===========================

%         % enforce measurement on exit waves:
%     %     phi( :, :, :, ss ) = enforce_2DTPAmeas( phi( :, :, :, ss ), meas_D( :, :, ss ), 1, sol );
%         phi( :, :, :, ss ) = enforce_2DTPAmeas( phi( :, :, :, ss ), expt.meas.SI( ss ), 1, sol );


meas_Deq0 = ( meas_D( :, :, ss ) == 0 );

     
V = fft( fftshift( fft( fftshift( phi( :, :, :, ss ), 1 ), [], 1 ), 2 ), [], 2 ) / sqrt_rc;
tmp0 = sqrt( sum( abs( V ) .^ 2, 3 ));

% exit waves corresponding to the different probe modes:
for pp = 1 : Nscpm             % size( phi, 3 )
  
    %==============
    
%     tmp1 = meas.D .* (( V( :, :, pp ) ./ ( 1e-7 + tmp0 )) .* meas.Dneq0 ) + V( :, :, pp ) .* meas.Deq0;
    tmp1 = meas_D( :, :, ss ) .* ( V( :, :, pp ) ./ ( 1e-7 + tmp0 )) + V( :, :, pp ) .* meas_Deq0;
    
    %==============
    
    phi( :, :, pp, ss ) = fftshift( ifft2( tmp1 )) * sqrt_rc;
    
%     V( :, :, pp ) = fftshift( ifft2( measLPF .* meas.D .* ( V( :, :, pp ) ./ ( 1e-7 + tmp0 )))) * sol.sz.sqrt_rc;    

    %==============
    
end  










    end

end
























function [ sol, expt ] = set_defaults( sol, expt )

    warning('off','MATLAB:prnRenderer:opengl');
    
    rng( 'shuffle' )
 
    %================================================================
    %--------- Gaussian LPF for reconstruction resolution -----------
    %================================================================
        
    mu    = 0.5 * sol.sz.sz + 1;
    stdev = 0.80 * sol.sz.sz;
    tmp0  = make_2Dgaussian( sol.sz.sz, mu, stdev );
    
    sol.measLPF = 1 + 0 * fftshift( tmp0 );

    %============================================
    %---------- Fixed sample support ------------
    %============================================
    
%     sol.sample.support = 0 + 1 * make_rectangle( sol.sample.sz.sz, [ 0.43 * sol.sample.sz.r, 0.43 * sol.sample.sz.c ]);
%     % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sample.sz.sz, [ 0.75 * sol.sample.sz.c, 0.75 * sol.sample.sz.r ]);
% 
%     [ sol.sample.support ] = 1 + 0 * lpf_gauss( sol.sample.support, 222.01 * sol.sample.sz.sz );
 
    %{
    
    figure; imagesc( sol.sample.support .* angle( sol.sample.TF ))
    figure; imagesc( sol.sample.support )
    
    %}
    

    %======================================================================
    %---------- Sample magnitude/phase scaling + modifications ------------
    %======================================================================
    
    %{
    
    abs_TF = abs(  sol.sample.TF );
    abs_TF = abs_TF / max( abs_TF( : ));
    
    phs_TF = angle(  sol.sample.TF );
    phs_TF = phs_TF - min( phs_TF( : ));
    phs_TF = phs_TF / max( phs_TF( : ));
    
    as = 0.2;
    ps = 0.5;
    
    tmp0 = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
    
    sol.sample.TF = tmp0;

    %}



    %{

    abs_TF = abs( sol.sample.TF );
    phs_TF = angle( sol.sample.TF );

    tmp0 = abs_TF .* exp( 1i * phs_TF + 1i * 0.0 );

    figure; 
    
    ax1 = subplot(131); imagesc( abs( tmp0 ));   daspect([1 1 1]); colormap( ax1, expt.cm.blkgrn )
    ax2 = subplot(132); imagesc( angle( tmp0 )); daspect([1 1 1]); colormap( ax2, expt.cm.blj )
    subplot(133);       imagescHSV( tmp0 );      daspect([1 1 1]);

    sol.sample.TF = tmp0;

    %}

    %============================================
    %---------- Fixed probe support -------------
    %============================================

    sol.probe.support = 0 + 1 * make_rectangle( sol.sz.sz, [ 0.9 * sol.sz.r, 0.9 * sol.sz.c ]);
    % sol.probe.support = 0 + 1 * make_ellipsoid( sol.sz.sz, [ 0.75 * sol.sz.c, 0.75 * sol.sz.r ]);

    [ sol.probe.support ] = 1 + 0 * lpf_gauss( sol.probe.support, 0.03 * sol.sz.sz );

    %======================================================
    %------ When to start sample/probe/spos updates -------
    %======================================================
    
    sol.it.spos_start  = 1e99;
    sol.it.spos_update = 1e99;
    
    sol.spos.suofai  = logical( 0 );                 % same update order for all iterations
    sol.spos.rand_ss = 1/2;                         % 

    sol.it.sample_start        = 0;
    sol.it.sample_update       = 1;
    sol.it.sample_sparsity     = 1;
    sol.it.sample_mag_phs_ineq = 1;
    
    sol.it.probe_start   = 0;
    sol.it.probe_update  = 1;
    sol.it.probe_orthog  = 1;
    sol.it.probe_scaling = 1;
    sol.it.probe_maxvals = 1;
    sol.it.probe_support = 1;

    sol.it.print_img_results = 10;
    sol.it.collect_metrics   = 50;

    %============================================
    %-------- PR projection algo params ---------
    %============================================
     
    sol.RAAR.beta = 0.5;

    %============================================
    %----- Sparsity constraints for sample ------
    %============================================
    
    sol.sparse.s2DFDxy = setup_edgedetect_forwarddiff2Dxy( [], sol.sample.sz.sz );
    
    sol.sparse.type = 'sparse2DFDxy';
    
    sol.sparse.type_sparse2DFDxy = { 'aniso_abs', ...                  % 1
                                     'aniso_phs', ...                  % 2
                                     'iso_abs', ...                    % 3
                                     'iso_phs', ...                    % 4
                                     'iso_abs_iso_phs', ...            % 5
                                     'aniso_abs_aniso_phs', ...        % 6
                                     'aniso_re_aniso_im', ...          % 7
                                     'iso_re_iso_im' };                % 8

    %============================================
    %-------------- Sample scaling --------------
    %============================================
    
    sol.sample.phsL = 0;
    sol.sample.phsH = +1.9 * pi;

    sol.sample.absL = 0.0;
    sol.sample.absH = 1.0;
    
    %============================================
    %-------------- Probe scaling ---------------
    %============================================
    
    %%{

    sol.probe.P = orthog_modes_eigendecomp( sol.probe.P );
    
    %========================

    sol.probe.scpm.mean_meas_fro2TOT = mean( expt.meas.SI_sumD2 );
%     sol.probe.scpm.fro2TOT = 4.20 * sol.probe.scpm.mean_meas_fro2TOT;
%     sol.probe.scpm.fro2TOT = 7e4 ^ 2;
    sol.probe.scpm.fro2TOT = 1400 ^ 2; 

    %========================
    
    sol.probe.scpm.max = 10 * [ 10, 20, 50 ];

    %========================
    
    sol.probe.scpm.occ = [ ];
%     sol.probe.scpm.occ = [ 0.05, 0.15, 0.80 ];
%     sol.probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
%     sol.probe.scpm.occ = [ 0.03, 0.07, 0.90 ];
    sol.probe.scpm.occ = sort( sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 ));     % make sure the mode occupancy adds up to 1.0
    
    sol.probe.P = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
    
    %}
    
    %============================================
    %---------- Scan positions init -------------
    %============================================
    
    sol.spos.update.dpx        = 1/4;  % for computing finite differences wrt scan position
    sol.spos.update.maxpx      = 4 * sol.spos.update.dpx;
    sol.spos.update.linesearch = transpose( 0.0 : sol.spos.update.dpx : sol.spos.update.maxpx );
    
%     sol.spos.update.shifttype = 'px';
    sol.spos.update.shifttype = 'subpx';

%     sol.spos.update.grad = 'all';
    sol.spos.update.grad = 'indiv';

    % for this data set, we use 0.5 um steps...only allow us to go so far from initial spos we start at
    sol.spos.update.maxcorrectr = 1.00e-6 / expt.csys.z2.dLv;
    sol.spos.update.maxcorrectc = 1.00e-6 / expt.csys.z2.dLh;
    
    sol.spos.indxsubset = expt.spos.indxsubset;                                                 % the scan positions we update the exitwaves/sample/probes over
    sol.spos.shifttype = 'px';                                                                  % When extracting part of the sample in the current scan position, use this type of pixel shifting  

    %================================================================
    %------ Center the probe function using center of mass ----------
    %================================================================
    
    %{

    figure; imagesc(abs(sol.probe.P( :, :, end )))

    [ com ] = centerOfMass( abs( sol.probe.P( :, :, end )));

    for pp = 1 : sol.probe.scpm.N

        sol.probe.P( :, :, pp ) = circ shift( sol.probe.P( :, :, pp ), -1 * round( com - 0.5 * sol.sz.sz - 1 ));

    end

    figure; imagesc(abs(sol.probe.P( :, :, end )))


    [~, Ic ] = max( sum( abs(sol.probe.P( :, :, end )), 1 ));
    [~, Ir ] = max( sum( abs(sol.probe.P( :, :, end )), 2 ));

    sol.probe.P( :, :, pp ) = nocircshift2D( sol.probe.P( :, :, pp ), -1 * round( [ Ir, Ic ] - 0.5 * sol.sz.sz - 1 ));

    figure; imagesc(abs(sol.probe.P( :, :, end )))

    %}


    %============================================
    %--------- Compute initial exit waves  ------
    %============================================

    
    
%     [ sol.phi ]       = ERupdate_2DTPAexwv_spos_meas( sol, expt );  
%     [ sol.sample.TF ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol );
%     [ sol.probe.P ]   = DMupdate_probemodes( sol.sample.TF, sol.phi, sol ); 
%     
%     sol.sample.TF = modulus_limits_project( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     % sol.sample.TF = modulus_limits_scale( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     sol.sample.TF = phase_limits_project( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );   
            


    sol.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N, sol.spos.N ], 'single' );

    for ss = 1 : sol.spos.N

        rs = +sol.spos.rs( ss, : );
    %     rs = -sol.spos.rs( ss, : );

        [ sol.phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +1 * rs, sol.spos.shifttype );

    end

    
    
    

end











































% 
% 
% 
% 
% 
% 
% 
% tic
% [ sample ] = DMupdate_sample( expt.probe.P, expt.phi, expt );
% toc
% 
% expt.spos.shifttype = 'px';
% tic
% [ TF ] = ePIEupdate_sample( expt.probe.P, expt.sample.TF, expt.phi, expt, expt );
% toc
% 
% 
% 
% 
% 
% 
% 
% expt.spos.updateorder = randperm( length( expt.spos.indxsubset ), round( 1.0 * length( expt.spos.indxsubset )));
%     
% 
% 
% 
% 
% 
% 
% 
% 







