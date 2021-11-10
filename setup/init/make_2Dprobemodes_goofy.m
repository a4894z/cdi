function [ probe, csys, sz ] = make_2Dprobemodes_goofy( expt )

%
% This is mostly just for quick and dirty algo testing purposes.
%
% Just define whatever goofy shapes you want to use for the probe function
% at the sample plane z2, and use a simple fft2 to numerically propagate
% the wave from the sample plane z2 to the measurement plane z3
%
%
%   |           |
%   |           |
%   |           |
%   |           |
%   |           |
%   |           |
%   |           |
%   z2          z3
%   sample      measurement
%   plane       plane
%

%====================================================================================================================================================
% -------------------------------------------------------- load previously defined SCPMs ------------------------------------------------------------
%====================================================================================================================================================

Z = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA_0.mat', 'expt' );

csys  = Z.expt.csys;
sz    = Z.expt.sz;
probe = Z.expt.probe;

probe.scpm.fro2TOT = 1 * probe.scpm.fro2TOT; 

% orthogonalize the probe modes:
[ probe.phi ] = orthog_modes_eigendecomp( probe.phi );

% ensure modes have the desired occupancy:
[ probe.phi, probe.scpm.fro2TOT, probe.scpm.occ ] = enforce_scpm_fro2TOT_photonocc( probe.phi, probe.scpm.fro2TOT, probe.scpm.occ );

return

%====================================================================================================================================================
% ---------------------------------------------------------------- create SCPMs  --------------------------------------------------------------------
%====================================================================================================================================================

%===========================================
% create coordinate system and problem sizes
%===========================================

csys.z23 = 5.00e-0;     % sample to detector ( in meters )

% pixel size of the 2D measurment device at z3 ( in meters )
csys.z3.dLy = 55e-6 / 1;    
csys.z3.dLx = 55e-6 / 1;

%========

sz.r = round( 256 + 0 ); 
sz.c = round( 128 + 0 );
sz.sz = [ sz.r, sz.c ]; 
sz.rc = sz.r * sz.c;
sz.sqrt_rc = sqrt( sz.rc );

%========

% fov of sample plane z2 using measurement plane z3 parameters
csys.z2.Ly = single( expt.lambda * csys.z23 / csys.z3.dLy ); 
csys.z2.Lx = single( expt.lambda * csys.z23 / csys.z3.dLx ); 
csys.z2.dLy = csys.z2.Ly / sz.r;
csys.z2.dLx = csys.z2.Lx / sz.c;

%===============================================================================
% fro^2 norm of all probe modes together, i.e. |P_1|^2 + |P_2|^2 + ... + |P_N|^2
%===============================================================================

probe.scpm.fro2TOT = 3000 ^ 2;  

%==========================================================
% (s)patially (c)oherent (p)robe (m)ode = SCPM occupancies:
%==========================================================

% probe.scpm.occ = [ 1.0 ];

% probe.scpm.occ = [ 0.2, 0.8 ];

% probe.scpm.occ = [ 0.32, 0.33, 0.35 ];
% probe.scpm.occ = [ 0.05, 0.20, 0.75 ];
% probe.scpm.occ = [ 0.05, 0.10, 0.85 ];
% probe.scpm.occ = [ 0.03, 0.07, 0.90 ];
% probe.scpm.occ = [ 0.02, 0.03, 0.95 ];

% probe.scpm.occ = [ 0.005, 0.015, 0.02, 0.03, 0.05, 0.08, 0.80 ];
% probe.scpm.occ = [ 0.006, 0.0085, 0.01, 0.0255, 0.03, 0.05, 0.87 ];

% probe.scpm.occ = [ 0.015, 0.025, 0.06, 0.12, 0.78 ];
% probe.scpm.occ = [ 0.01, 0.025, 0.06, 0.12, 0.78 ];
% probe.scpm.occ = [ 0.02, 0.03, 0.1, 0.15, 0.7 ];
% probe.scpm.occ = [ 0.2, 0.2, 0.2, 0.2, 0.2 ];

% probe.scpm.occ = random('unif', 0, 1, length( paths.probe ));

probe.scpm.occ = exp( -4 * linspace( 1, 0, 5 ));

%========

% make sure the mode occupancy adds up to 1.0:
probe.scpm.occ = sort( probe.scpm.occ / norm( probe.scpm.occ, 1 ));

% get total number of scpm used to define simulated x-ray beam:
probe.scpm.N = length( probe.scpm.occ );

%======================================================
% create gaussian like probe modes (not orthogonal yet)
%======================================================

probep2 = zeros( [ sz.sz, probe.scpm.N ], 'single' );

for pp = 1 : probe.scpm.N
    
%     tmpz0 = make_rectangle( sz.sz, ( 0.1 * ( 2 * rand  - 1 ) + 0.4 ) * [ min( sz.sz ), min( sz.sz ) ] );
%     tmpz0 = make_rectangle( sz.sz, ( 0.2 * ( 2 * rand  - 1 ) + 0.3 ) * sz.sz );
%     tmpz0 = make_rectangle( sz.sz, ( 0.1 * ( 2 * rand  - 1 ) + 0.3333 ) * sz.sz );

%     tmpz0 = make_2Dellipsoid( sz.sz, ( 0.3 * ( 2 * rand  - 1 ) + 0.3 ) * sz.sz );
%     tmpz0 = make_2Dellipsoid( sz.sz, 0.6 * sz.sz );

    tmpz0 = make_2Dgaussian( sz.sz, ( 0.5 *  sz.sz + 1 ), [ ( 0.01 * ( 2 * rand  - 1 ) + 0.06 ) * sz.sz( 1 ), ...
                                                            ( 0.01 * ( 2 * rand  - 1 ) + 0.06 ) * sz.sz( 2 ) ]);


%     tmpz0 = image2complex( paths.probe{ pp } );
%     tmpz0 = padarray( tmpz0, [ 100, 100 ] );
%     tmpz0 = imresize( tmpz0, sz.sz, 'nearest' );  % box, triangle, bilinear, nearest

    tmpz1 = lpf_gauss( tmpz0, 0.005e0 * sz.sz ); 
    tmpz2 = lpf_gauss( tmpz0, 0.010e0 * sz.sz );
    
    probep2( :, :, pp ) = tmpz2 .* exp( 1i * 2 * pi * 3.0 * tmpz1 );
    
    tmp0 = max( max( abs( probep2( :, :, pp ))));
    probe.phi( :, :, pp ) = ( abs( probep2( :, :, pp )) > 1e-3 * tmp0 ) .*  probep2( :, :, pp );
    
%     probe.phi( :, :, pp ) = probe.phi( :, :, pp ) .* ( probe.phi( :, :, pp ) > 1e-5 );
    
end






%====================================================================================================================================================

% ensure modes have the desired occupancy:
[ probe.phi, probe.scpm.fro2TOT, probe.scpm.occ ] = enforce_scpm_fro2TOT_photonocc( probe.phi, probe.scpm.fro2TOT, probe.scpm.occ );

%========

% for pp = 1 : probe.scpm.N
%    
%     figure; 
%     subplot(121); imagesc( log10(1 + 1e0 * abs( probe.phi( :, :, pp ) ))); daspect([1 1 1])
%     subplot(122); imagescHSV( log10(1 + 1e0 * abs( probe.phi( :, :, pp ) )) .* exp(1i * angle(probe.phi( :, :, pp )))); daspect([1 1 1])
% end
% 
% close all;

%====================================================================================================================================================

% orthogonalize the probe modes:
[ probe.phi ] = orthog_modes_eigendecomp( probep2 );

% ensure modes have the desired occupancy:
[ probe.phi, probe.scpm.fro2TOT, probe.scpm.occ ] = enforce_scpm_fro2TOT_photonocc( probe.phi, probe.scpm.fro2TOT, probe.scpm.occ );


for pp = 1 : probe.scpm.N
    
    figure; 
    subplot(121); imagescHSV( log10(1 + 10^0 * abs( probe.phi( :, :, pp ) )) .* exp(1i*angle(probe.phi( :, :, pp ))) ); daspect([1 1 1])
    subplot(122); imagesc( log10(1 + 10^0 * abs( probe.phi( :, :, pp ) ))); daspect([1 1 1]); colormap( expt.cm.blj ); colorbar
    
end

figure; imagesc( sum( abs( probe.phi ) .^ 2, 3 )); daspect([1 1 1]); colormap( expt.cm.blj ); colorbar

close all;

%====================================================================================================================================================

%{

%====================
% check orthogonality
%====================

tmp1 = reshape( probe.phi, [ sz.r * sz.c, probe.scpm.N ] );

for ii = 1 : probe.scpm.N

    tmp0 = probe.phi( :, :, ii );
    tmp0 = tmp0( : );
    
    tmp2( ii, : ) = tmp0' * tmp1;

end

figure; imagesc(sqrt(abs(tmp2)))

%}


