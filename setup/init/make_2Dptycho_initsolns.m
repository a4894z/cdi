function [ sol ] = make_2Dptycho_initsolns( expt )

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

sol.sample.T = lpf_gauss( sol.sample.T, [ 0.1 * expt.sample.sz.r, 0.1 * expt.sample.sz.r ] );

sol.sample.T = modulus_limits_scale( sol.sample.T, [ 0.01, 0.999 ] );
sol.sample.T = phase_limits_scale( sol.sample.T, [ 0 * pi, 0.7 * pi ] );

sol.sample.sz = expt.sample.sz;

%==============
% 2D SCPM probe 
%==============

% some goofy shapes for the probe modes

sol.probe.scpm = expt.probe.scpm;

% sol.probe.scpm.occ = [ 1 ];
% sol.probe.scpm.occ = [ 0.05, 0.1, 0.2, 0.4, 0.9 ];
% sol.probe.scpm.occ = [ 0.2, 0.8 ];
% sol.probe.scpm.occ = [ 0.01, 0.04, 0.95 ];
% sol.probe.scpm.occ = [ 0.32, 0.33, 0.35 ];
% sol.probe.scpm.occ = [ 0.2, 0.2, 0.2, 0.2, 0.2 ];
sol.probe.scpm.occ = exp( -7 * linspace( 1, 0, 8 ));

sol.probe.scpm.occ = sol.probe.scpm.occ / norm( sol.probe.scpm.occ, 1 );

sol.probe.scpm.N = length( sol.probe.scpm.occ );

guess_probe = zeros( [ sol.sz.sz, sol.probe.scpm.N ], 'single' );

for pp = 1 : sol.probe.scpm.N
    
    mmu     = [ 0.5 * expt.sz.r + 1, 0.5 * expt.sz.c + 1 ];
%     sstev   = [ 0.20 * expt.sz.r, 0.20 * expt.sz.c ];
    sstev   = [ 0.02 * ( 2 * rand + 2 ) * expt.sz.r, 0.02 * ( 2 * rand + 1 ) * expt.sz.c ];
    
    tmp0 = make_2Dgaussian( expt.sz.sz, mmu, sstev );
    
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

sol.spos.rs = expt.spos.rs;

sol.spos.N = length( sol.spos.rs );

% sol.spos.indx = 1 : expt.spos.N;
sol.spos.indx = expt.spos.indx;

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

