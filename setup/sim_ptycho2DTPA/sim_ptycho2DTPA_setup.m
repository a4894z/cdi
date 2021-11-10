function [ sol, expt ] = sim_ptycho2DTPA_setup( expt )

%
% TODO: 
% - make probe and sample at higher resolution, e.g. 1 nm then let measurement problem
% size determine final problem size (i.e. down interpolate the calculated measurement based
% on some desired detector pixel size.
%
% -
%
% -
% setup_simPRexpt_forwardmodel2DTPA
% setup_PRexpt_forwardmodel2DTPA
% setup_forwardmodel_2DTPA
%

%{

clear; close all; 
restoredefaultpath
% expt.path.code = '~/Documents/Science/Matlab/Code/cdi/';   
% cd( expt.path.code );
addpath(genpath(pwd));





clear; close all; [ sol, expt ] = sim_ptycho2DTPA_setup






[ V ] = make_gaussian( [512,512], [257, 257], 0.2 * [512,512] ); 
V = V * 1e3 / norm( V, 'fro' );
figure; imagesc( V )
figure; imagesc( max( abs( V(:) ).^2 ) - abs( V ).^2 )



1e-6 m probe
75% overlap --> 0.25e-6 shifts
0.25e-6 *


%}

rng( 'shuffle' )

%==================================================================================================
%---------------------------------------------- BEGIN ---------------------------------------------
%==================================================================================================

fprintf('\n========================================================================================================'); 
fprintf('\nCreating Simulated Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed...\n'); 
fprintf('========================================================================================================\n'); 
    
% expt = struct;
% sol  = struct;

%======================================
% set up paths for phase retrieval code
%======================================

% expt.paths.code     = '~/Documents/Science/Matlab/Code/cdi/';   
% expt.paths.rimgdata = '/media/ash/Saltwater/Data/simulated/';   

% expt.paths.code     = '~/Documents/Science/Matlab/Code/cdi/';   
% expt.paths.rimgdata = '/home/ash/Documents/Science/Matlab/Data/simulated/img/';   


expt.paths.code     = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/';   
expt.paths.rimgdata = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/simulated/img/';   

%================
% load misc stuff 
%================

% color maps I like:

tmp1 = load( [ expt.paths.code, 'misc/colormaps/colormap_blue_light_jet.mat' ]);
expt.cm.blj = tmp1.blue_light_jet;
clear( 'tmp1' );

tmp1 = load( [ expt.paths.code, 'misc/colormaps/colormap_black_green.mat' ]);
expt.cm.blkgrn = tmp1.blkgrn;
clear( 'tmp1' );

tmp1 = load( [ expt.paths.code, '/misc/colormaps/colormap_hsvD.mat' ]);
expt.cm.hsvD = tmp1.hsvD;
clear( 'tmp1' );

%====================================================
% define the various relevant experimental parameters 
%====================================================

expt.energy = 10.000;                           % in keV
expt.lambda = ( 12.4 / expt.energy ) * 1e-10;   % wavelength ( in meters )

%==================================================
% load/create/define spatially coherent probe modes 
%==================================================

% % load from images, past results, from experiment, etc:
% [ expt.probe, expt.csys, expt.sz ] = make_2Dprobemodes_loadfromfile( expt );

%========

% using whatever goofy shapes you like:
[ expt.probe, expt.csys, expt.sz ] = make_2Dprobemodes_goofy( expt );

%========

% % using a beam defining aperture ( e.g. pinhole ):
% [ expt.probe, expt.csys, expt.sz ] = make_2Dprobemodes_UFB( expt );

%========

% % using a focusing optic:
% [ expt.probe, expt.csys, expt.sz ] = make_2Dprobemodes_FB( expt );

%========

% clean up:
expt      = orderfields( expt );
expt.csys = orderfields( expt.csys );

%==================================
% load/create/define scan positions
%==================================

% DESCRIBE WHAT WE'RE DOING HERE
% DESCRIBE WHAT WE'RE DOING HERE
% DESCRIBE WHAT WE'RE DOING HERE
% DESCRIBE WHAT WE'RE DOING HERE

% percent overlap is: ( 1 - step size / expt.probe diameter )

% define downstream optical axis as +y, then +x ( horizontal ) is inboard 
% (towards storage ring as viewed from beamline ), and +z (vertical) is up.

% !!!!!!!!!!!!!!! spos should be [ index, row pos, col pos ]    !!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!! currently spos is [ row pos, col pos ]        !!!!!!!!!!!!!!!

%=====================================
% scan positions on a rectangular grid
%=====================================

% [ expt.spos ] = make_2Dptycho_spos_rectangular( expt );

%=====================================
% scan positions on concentric circles
%=====================================

[ expt.spos ] = make_2Dptycho_spos_concentric_circles( expt );

%=================
% plotting of spos
%=================

% tmp0 = randperm( expt.spos.N, round( 0.05 * expt.spos.N) );
% 
%     % plot scan positions
% 
%     figure; 
%     plot_2Dscan_positions( expt.spos.rs, [], expt.spos.rs( tmp0, : ), [] )
% %     plot_2Dscan_positions( expt.spos.yv_xh_frameL, [], expt.spos.rs( 1 : 1 : end, : ), expt.spos.indxsubset( 1 : 1 : end ))
%     set( gca, 'xdir', 'reverse' )
%     set( gca, 'ydir', 'normal' )
%     daspect([1 1 1])
% %     xlabel('xh, lab frame'); ylabel('yv, lab frame');
%     %title('positions for scanning the probe on the sample')
% 
%    
%     
%     
% 
%     % distance between adjacent scan positions (in spos order)         
% 
%     tmp = zeros( expt.spos.N - 1, 1, 'single' );
% 
%     for ss = 1 : expt.spos.N - 1
% 
%         tmp( ss ) = sqrt(( expt.spos.rs( ss, 1 ) - expt.spos.rs( ss + 1, 1 )) .^ 2 + ( expt.spos.rs( ss, 2 ) - expt.spos.rs( ss + 1, 2 )) .^ 2);
% 
%     end
% 
%     figure
%     plot( tmp, 'Linewidth', 1 ); 
%     xlabel('scan position index'); ylabel('distance between adjacent scan positions')
%     grid on


%===============================
% introduce scan position errors
%===============================

% ???

%=========
% clean up
%=========

expt      = orderfields( expt );
expt.spos = orderfields( expt.spos );
% clearvars -except expt

%=======================================
% load/create/define sample transmission 
%=======================================

[ expt.sample ] = make_2Dsample( expt );

%========

% %load ~/Documents/MATLAB/Code/misc/sparsity_paper_code_15Aug2016/simulated_experiment/test_problems/optexp_test_probs/512_16e.mat
% load ~/Documents/MATLAB/Code/misc/sparsity_paper_code_15Aug2016/simulated_experiment/test_problems/optexp_test_probs/512_512e.mat
% expt.sample.T = smpl_trans;
% 
% expt.sample.T = modulus_limits_scale( expt.sample.T, [ 0.5, 0.9 ] );
% expt.sample.T = phase_limits_scale( expt.sample.T, 2 * pi * [ -4, 4 ] );
% 
% expt.sample.sz.r = single( size( expt.sample.T, 1 ));
% expt.sample.sz.c = single( size( expt.sample.T, 2 ));
% expt.sample.sz.rc = expt.sample.sz.r * expt.sample.sz.c;
% expt.sample.sz.sqrt_rc = sqrt( expt.sample.sz.rc );
% expt.sample.sz.sz = [ expt.sample.sz.r, expt.sample.sz.c ]; 
% 
% % view slices for when computing the exit wave views at a particular scan position
% expt.sample.vs.r = round( ( 0.5 * ( expt.sample.sz.r - expt.sz.r ) + 1 ) : ( 0.5 * ( expt.sample.sz.r + expt.sz.r )));
% expt.sample.vs.c = round( ( 0.5 * ( expt.sample.sz.c - expt.sz.c ) + 1 ) : ( 0.5 * ( expt.sample.sz.c + expt.sz.c )));

%========

% clean up 
expt = orderfields( expt );
expt.sample = orderfields( expt.sample );
clearvars -except expt sol

%========

% does_the_sample_array_fit_all_spos( expt )
% 
% close all;

%================================
% define diffraction measurements
%================================

% generate measurements according to which type of wavefield propagation 
% is used to get from the sample plane z2 to the measurement plane z3.

%========

% assuming an unfocused beam is used (e.g. a pinhole ) or sample is at the focal plane:
[ expt.meas, expt.csys.z3 ] = make_measurements_2DTPA_UFB( expt );

expt.csys = orderfields( expt.csys );
expt.meas = orderfields( expt.meas );

%========

% % assuming sample is at a sample plane defocused ( not at focal plane z1 ):
% [ expt.meas ] = make_measurements_2DTPA_FB( expt );

%======================================================================
% initial definitions for the parameters we wish to solve for solve for  
%======================================================================

[ sol ] = make_2Dptycho_initsolns( expt );

% counters keeping track of the total number of iterations we've 
% run, defined as an update over all or a subset of scan locations. 
sol.it.mtot = 1;
sol.it.metr = 1;
sol.it.exwv = 1;
sol.it.epoch = 1;

%==================================================================================================
%--------------------------------------- final cleaning up ----------------------------------------
%==================================================================================================

clearvars -except expt sol
expt = orderfields( expt );
sol  = orderfields( sol );

%==================================================================================================
%-------------------------- create save file from quantities defined ------------------------------
%==================================================================================================

% expt.paths.rsdata = [ expt.paths.code, '/run/misctesting2DTPA/', 'simPRexpt_2DTPA.mat' ]; % saved data path
expt.paths.rsdata = [ expt.paths.code, 'sim_ptycho2DTPA.mat' ]; % saved data path

fprintf('saving simulated expt ...'); 

save( expt.paths.rsdata, '*'); 
copyfile( expt.paths.rsdata, [ expt.paths.rsdata( 1 : end - 4 ), '_0' , '.mat'] );

fprintf(' done!\n\n');

%==================================================================================================
%---------------------------------------------- END -----------------------------------------------
%==================================================================================================

end

function does_the_sample_array_fit_all_spos( sol )


    sol.spos.shifttype = 'px';
    % expt.spos.shifttype = 'subpx';


    pm = size( sol.probe.P, 3 );


    tmp0 = 1 + 0 * single( abs(sol.probe.P( :, :, pm )) > 1e-9 );
    
%     max_P = max(max( abs(sol.probe.P( :, :, pm ))));
%     tmp0 = single( abs(sol.probe.P( :, :, pm )) > 1e-4 * max_P );
    % tmp9 = padarray( tmp0, double( 0.5 * [ sol.sample.sz.r - expt.sz.r, sol.sample.sz.c - expt.sz.c ] ));
    tmp9 = zeropadarray( tmp0, double( 0.5 * [ sol.sample.sz.r - sol.sz.r, sol.sample.sz.c - sol.sz.c ] ));

    % % just used as a visual aide to see where ~ center of probe is
    % tmp9( size( tmp9, 1 ) / 2 + 1, size( tmp9, 2 ) / 2 + 1 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 0, size( tmp9, 2 ) / 2 + 1 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 0, size( tmp9, 2 ) / 2 + 0 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 1, size( tmp9, 2 ) / 2 + 0 ) = 10;  

    % tmp0 = ( abs(sol.sample.T) > 0.01 * max(abs(sol.sample.T(:))) );
    tmp55 = 0;

    for ii = 1 : sol.spos.N

    %     rj = round( -spos.rs( sio( ii ), : ));
        rj = round( +sol.spos.rs( ii, : ));

        tmp2 = circshift( round( tmp9 ), -rj );

        [ phi, ~ ] = enforce_2DTPAsposview( sol.probe.P( :, :, pm ), ( 0 + 1 * sol.sample.T ), sol.sample.vs, +rj, sol.spos.shifttype );

        tmp55 = tmp55 + tmp2 .* sol.sample.T;


        %{

        figure(666);
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1] );

    %     tmp0 = expt.meas.SI( ii ).Dneq0 .* expt.meas.SI( ii ).D + expt.meas.SI( ii ).Deq0 .* fft2( phi ) / sol.sz.sqrt_rc;
        tmp0 = fft2( phi ) / sol.sz.sqrt_rc;

    % p(1)=subplot(3,4,[1 2 5 6]);
    % p(2)=subplot(3,4,[3 4 7 8]);
    % p(3)=subplot(3,4,[9 10]);
    % p(4)=subplot(3,4,[11 12])

        subplot( 2, 3, 1 ); imagesc( log10(1 + abs( fftshift( tmp0 )))); daspect([1 1 1]); colormap jet; colorbar
        subplot( 2, 3, 4 ); imagesc( log10(1 + abs( fftshift( expt.meas.SI( ii ).D )))); daspect([1 1 1]); colormap jet; colorbar
        %subplot(122); imagesc( ( tmp1 > 0 ) .* abs( 1 + 0 * sol.sample.T ) ); daspect([1 1 1]); colormap jet
        %subplot(132); imagesc( abs(phi) > 0 ); daspect([1 1 1]); colormap jet
        subplot( 2, 3, [ 2, 5 ] ); imagescHSV( phi + 0 * tmp0 ); daspect([1 1 1]); %colormap jet
        subplot( 2, 3, [ 3, 6 ] ); imagescHSV( tmp2 .* sol.sample.T ); daspect([1 1 1]); %colormap jet
    %     tic
    %     export_fig( num2str( ii, 'test_spos_%d' ), '-r90.0' )
    %     toc
        pause( 0.01 )

        %}


    end

    figure; imagesc(abs(tmp55)>0); colormap gray; daspect([1 1 1])

    % figure; imagescHSV( ( tmp1 > 0 ) .* ( 1 + 0 * sol.sample.T )); daspect([1 1 1]); %colormap gray


%     close all;

end










%{

clear; close all; 
restoredefaultpath
expt.paths.code = '~/Documents/Science/Matlab/Code/cdi/';   
cd( expt.paths.code );
addpath( genpath( pwd ));
clearvars -except expt



clear; close all; ptycho2DTPA_setup;








tic
for i = 1:4
    v = zeros(1, 10);
    for j = 1:10
        v(j) = i + j;
    end
    disp(v)
    A(i, :) = v;
end
toc

tic
parfor i = 1:4
    v = zeros(1, 10);
    for j = 1:10
        v(j) = i + j;
    end
    disp(v)
    A(i, :) = v;
end
toc














% spos = [ [ 0, 0 ]; [ 0, 1 ]; [ 1, 0 ]; [ 1, 1 ]; [ 2, 0 ]; [ 2, 1 ] ];
% spr = repmat( 0:szr-1, szc, 1 ); spr = spr(:)';
% spc = repmat( 0:szc-1, 1, szr ); Tc = Tc(:)';








clear;

Nr = 6;
Nc = 4;
Np = 10;

Av = 1 * ( 1 : ( Nr * Nc * Np ));
Am = reshape( Av, [ Nr, Nc, Np ] );


% these must always be even
szr = 2;
szc = 3;

sz = [ szr, szc ];

vsrb = round( ( 0.5 * Nr + 1 - 0.5 * szr ) : ( 0.5 * Nr + 0.5 * szr ));
vscb = round( ( 0.5 * Nc + 1 - 0.5 * szc ) : ( 0.5 * Nc + 0.5 * szc ));





vsr = repmat( 1 : szr, szc, 1 ); 
vsc = repmat( 1 : szc, 1, szr ); 
vsr = vsr(:)';
vsc = vsc(:)';







spos( :, 1 ) = [ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4 ];
spos( :, 2 ) = [ 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 ];



s = 2;
spr = ( 1 : sz( 1 )) + spos( s, 1);
spc = ( 1 : sz( 2 )) + spos( s, 2);

Am( :, :, 1 )
Am( spr, spc, 1 )



return












Ir = [ vsr + 0, vsr + 0, vsr + 1, vsr + 1, vsr + 2, vsr + 2, vsr + 3, vsr + 3, vsr + 4, vsr + 4 ];
Ic = [ vsc + 0, vsc + 1, vsc + 0, vsc + 1, vsc + 0, vsc + 1, vsc + 0, vsc + 1, vsc + 0, vsc + 1 ];
Ip = repmat( 1:Np, szr * szc, 1 ); Ip = Ip( : )';


ind = sub2ind( [ Nr, Nc, Np ], Ir, Ic, Ip )





tic, tt = Am( ind ), toc
ff = Am(:);
tic, ff( ind ), toc










return




%}





