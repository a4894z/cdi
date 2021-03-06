function [ sol ] = make_2Dptycho_initsolns( expt )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD PREVIOUSLY DEFINED CANDIDATE SOLUTION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Z = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA.mat', 'sol' );
% 
% sol = Z.sol;
% 
% sol.spos.indxsubset = expt.spos.indxsubset;
% 
% %========
% 
% sol.probe.scpm.fro2TOT = expt.probe.scpm.fro2TOT;
% % sol.probe.scpm.fro2TOT = ( 0.2 * ( 2 * rand - 1 ) + 1.0 ) * expt.probe.scpm.fro2TOT;
% 
% % ensure modes have the desired occupancy:
% [ sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ ] = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );
% 
% %========
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE CANDIDATE SOLUTION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============================
% misc experimental parameters
%=============================

sol.csys               = expt.csys;
sol.lambda             = ( 0.0 * ( 2 * rand - 1 ) + 1.0 ) * expt.lambda;
sol.probe.scpm.fro2TOT = ( 0.2 * ( 2 * rand - 1 ) + 1.0 ) * expt.probe.scpm.fro2TOT;
sol.sz                 = expt.sz;

%==========
% 2D sample 
%==========

% blurred random complex numbers for the sample transfer function

sol.sample.T = rand( expt.sample.sz.r, expt.sample.sz.c ) .* exp( 1i * 2 * pi * rand( expt.sample.sz.r, expt.sample.sz.c ));

sol.sample.T = lpf_gauss( sol.sample.T, [ 0.02 * expt.sample.sz.r, 0.02 * expt.sample.sz.c ] );

sol.sample.T = modulus_limits_scale( sol.sample.T, [ 0.01, 0.999 ] );
sol.sample.T = phase_limits_scale( sol.sample.T, [ 0 * pi, 0.1 * pi ] );

sol.sample.sz = expt.sample.sz;

%==============
% 2D SCPM probe 
%==============

sol.probe.scpm.N = expt.probe.scpm.N + 2;

sol.probe.scpm.occ = exp( -4 * linspace( 1, 0, sol.probe.scpm.N ));
sol.probe.scpm.occ = sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 );

guess_probe = zeros( [ sol.sz.sz, sol.probe.scpm.N ], 'single' );

for pp = 1 : sol.probe.scpm.N
    
    sstev   = [ 0.02 * ( 2 * rand + 2 ) * expt.sz.r, ...
                0.02 * ( 2 * rand + 1 ) * expt.sz.c ];
    
    tmp0 = make_2Dgaussian( expt.sz.sz,             ...
                            0.5 * expt.sz.sz + 1,   ...
                            sstev );
    
    tmp0 = tmp0 / max( abs( tmp0(:) ));     
         
    tmp0 = tmp0 .* exp( 2 * pi * 1i * 0 * tmp0 );
    
    guess_probe( :, :, pp ) = tmp0 * sol.probe.scpm.fro2TOT / norm( tmp0, 'fro' );

end

% orthogonalize the probe modes:
[ sol.probe.phi ] = orthog_modes_eigendecomp( guess_probe );

% ensure modes have the desired occupancy:
[ sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ ] = enforce_scpm_fro2TOT_photonocc( sol.probe.phi, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ );

% for pp = 1 : sol.probe.scpm.N
%     
%     figure; 
%     subplot(221); imagescHSV( log10(1 + 10^0 * abs( sol.probe.phi( :, :, pp ) )) .* exp(1i*angle( sol.probe.phi( :, :, pp ))) ); daspect([1 1 1])
%     subplot(222); imagesc( log10(1 + 10^0 * abs( sol.probe.phi( :, :, pp ) ))); daspect([1 1 1]); colormap jet
%     subplot(223); imagescHSV( log10(1 + 10^0 * abs( expt.probe.phi( :, :, pp ) )) .* exp(1i*angle( expt.probe.phi( :, :, pp ))) ); daspect([1 1 1])
%     subplot(224); imagesc( log10(1 + 10^0 * abs( expt.probe.phi( :, :, pp ) ))); daspect([1 1 1]); colormap jet
%  
% end
%    
% close all;

%=========================
% 2D ptycho scan positions 
%=========================

sol.spos.shifttype = 'px';
% sol.spos.shifttype = 'subpx';  

sol.spos.rs = expt.spos.rs0;


% figure; 
% plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-600, 600])
% ylim([-600, 600])
% daspect([1 1 1])  
% grid on
% 
% 5;




sol.spos.N = length( sol.spos.rs );

sol.spos.indx       = expt.spos.indx;
sol.spos.indxsubset = sol.spos.indx;

%========

% view slices for when computing the exit wave views at a particular scan position
sol.sample.vs.r = round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r )));
sol.sample.vs.c = round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c )));

%===========================
% misc other stuff we'll use 
%===========================

% counters keeping track of the total number of iterations we've 
% run, defined as an update over all or a subset of scan locations. 
sol.it.mtot = 1;
sol.it.metr = 1;
sol.it.exwv = 1;

