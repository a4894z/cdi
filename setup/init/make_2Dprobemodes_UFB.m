function [ probe, csys, sz ] = make_2Dprobemodes_UFB( expt )

% make_2Dprobemodes_F
% make_2Dprobemodes_UF

%==================================================================================================
%
% We need to define a 2D pinhole transmission function at z0 plane,
% then propagate that wavefield to the sample location at z2. To do 
% that, we need to make sure we can use spectrum propagation:
%
%
%   |               |           |
%   |               |           |
%   |               |           |
%   |               |           |
%   |               |           |
%   |               |           |
%   |               |           |
%   z0              z2          z3
%   beam defining   sample      measurement
%   aperture        plane       plane    
%   plane

%==================================================================================================

csys.z02 = 1.00e-3;                % beam defining aperture plane z0 to sample plane z2 ( in meters )
%csys.z23 = 6.451612903226e-3;                 % sample plane z2 to measurement plane z3 ( in meters )
csys.z23 = 1e0;                 % sample plane z2 to measurement plane z3 ( in meters )
csys.z03 = csys.z02 + csys.z23;     % beam defining aperture plane z0 to measurement plane z3 ( in meters )

%================================================

paths.probe = cell( 0 );

% % paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'zp6b_512px.ppm' ];
% % % % paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'PETN2.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'pinholeA.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'zp6_512px.ppm' ];
% % paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'spiral_zp.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'pinholeC.ppm' ];

paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'pinholeD.ppm' ];
paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'PETN2.ppm' ];
paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'pinholeGa.ppm' ];

%================================================

% (s)patially (c)oherent (p)robe (m)ode occupancies:

% probe.scpm.occ = 1;

% % probe.scpm.occ = [ 0.3, 0.7 ];
% probe.scpm.occ = [ 0.5, 0.5 ];

probe.scpm.occ = [ 0.05, 0.10, 0.85 ];

% probe.scpm.occ = [ 0.25, 0.25, 0.5 ];
% probe.scpm.occ = [ 0.33, 0.33, 0.34 ];
% probe.scpm.occ = [ 0.1, 0.2, 0.7 ];+

% probe.scpm.occ = [ 0.25, 0.25, 0.25, 0.25 ];
% probe.scpm.occ = [ 0.05, 0.1, 0.2, 0.65 ];

% probe.scpm.occ = [ 0.1, 0.15, 0.2, 0.25, 0.3 ];
% probe.scpm.occ = [ 0.045, 0.065, 0.14, 0.25, 0.5 ];
% probe.scpm.occ = random('unif', 0, 1, length( paths.probe ));

% make sure the mode occupancy adds up to 1.0:
probe.scpm.occ = sort( probe.scpm.occ / norm( probe.scpm.occ, 1 ));

% get total number of scpm used to define simulated x-ray beam:
probe.scpm.N = length( probe.scpm.occ );

%================================================

% initial problem size for the beam defining aperture plane z0

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% this might be different than the final problem 
% size since we maybe will need to zero pad in order 
% to allow use of spectrum propagation from beam
% defining aperture plane z0 to sample plane z2 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sz.r = round( 200 ); 
sz.c = round( 200 );
sz.sz = [ sz.r, sz.c ]; 
sz.rc = sz.r * sz.c;
sz.sqrt_rc = sqrt( sz.rc );

%==================================================================================================

% pixel size of the beam defining aperture plane z0:
csys.z0.dr = 100.0e-9;
csys.z0.dc = 100.0e-9;

% focusing optic pupil function dimensions:
PFz0.Lr = 10e-6;  % vertical diameter/side length ( in meters )
PFz0.Lc = 10e-6;  % horizontal diameter/side length ( in meters )

PFz0.r = round( PFz0.Lr / csys.z0.dr );  % rectangular vertical length ( in pixels )
PFz0.c = round( PFz0.Lc / csys.z0.dc );  % rectangular horizontal length ( in pixels )

if PFz0.r > sz.r, error('Must increase the number of ROW pixels in the focusing optic pupil function plane z0.'); end
if PFz0.c > sz.c, error('Must increase the number of COLUMNS pixels in the focusing optic pupil function plane z0.'); end

clear( 'PFz0' );

%==================================================================================================

% create non-orthogonal modes at the beam defining aperture plane z0

probez0 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

abs_minmax_z0 = repmat( [ 0, 1 ], probe.scpm.N, 1);
phs_minmax_z0 = repmat( 1.3 * [ -pi, pi ], probe.scpm.N, 1);

for pp = 1 : probe.scpm.N
   
    tmpz0 = image2complex( paths.probe{ pp } );
    
%     tmpz0 = lpf_gauss( tmpz0, 0.05 * sz.sz ); 

    %tmpz0 = padarray( tmpz0, 2 * [ 16, 16 ] + 0 );
    tmpz0 = zeropadarray( tmpz0, 12 * [ 16, 16 ] + 0 );  
    
    %tmpz0 = imresize( tmpz0, sz.sz, 'bilinear' );        % box, triangle, bilinear, nearest
    tmpz0 = imresize2D( tmpz0, sz.sz );

    tmpz0 = modulus_limits_scale( tmpz0, abs_minmax_z0( pp, : ));
    tmpz0 = phase_limits_scale( tmpz0, phs_minmax_z0( pp, : ));
    
%     [ min( abs(tmpz0(:))), max( abs(tmpz0(:))) ]
%     [ min( angle(tmpz0(:))), max( angle(tmpz0(:))) ]

    probez0( :, :, pp ) = tmpz0;
    
end

% for pp = 1 : probe.scpm.N
%    
%     figure; 
%     subplot(121); imagescHSV( log10(1 + 1e0 * abs( probez0( :, :, pp ))) .* exp(1i * angle(probez0( :, :, pp )))); daspect([1 1 1])
%     subplot(122); imagesc( log10(1 + 1e0 * abs( probez0( :, :, pp ) ))); daspect([1 1 1]);  colormap(expt.cm.blkgrn);
% 
% end
% 
% close all;

%==================================================================================================
% zero pad the wavefield at the beam defining aperture plane z0 so we could use 
% spectrum propagation to numerically propagate the wavefield at the beam defining 
% aperture z0 to the sample plane z2

z0pad = round( [ 0, 0 ] + 0.5 * 56 );     % for sz.r, c = round( 200 )
%z0pad = round( [ 0, 0 ] + 0.5 * 28 );     % for sz.r, c = round( 100 )

sz.r = sz.r + ( 2 * z0pad( 1 ));
sz.c = sz.c + ( 2 * z0pad( 2 ));
sz.sz = [ sz.r, sz.c ]; 
sz.rc = sz.r * sz.c;
sz.sqrt_rc = sqrt( sz.rc );

csys.z0.Lr = csys.z0.dr * sz.r;
csys.z0.Lc = csys.z0.dc * sz.c;

%===============================

zcrit0.r = csys.z0.dr ^ 2 * ( 2 * z0pad( 1 )) / expt.lambda;
zcrit0.c = csys.z0.dc ^ 2 * ( 2 * z0pad( 2 )) / expt.lambda;

fprintf( '\nCheck to see if we can use Spectrum Propagation to get from z0 to z2:');
fprintf( [ '\n', num2str( 1e3 * [ zcrit0.r, zcrit0.c ], 'zcrit0 = ( %.3f, %.3f ) mm' ), ...
                 num2str( 1e3 * csys.z02, ', z02 = %.3f mm'), ...
                 num2str( 1e3 * csys.z23, ', z23 = %.3f mm'), ...
                 num2str( 1e3 * csys.z03, ', z03 = %.3f mm'), ' \n\n' ])

%===============================

probez0pad = zeros( sz.r, sz.c , 'single' );

% mu = [ sz.r / 2 + 1, sz.c / 2 + 1 ];
% stdev = [  0.05 * sz.r, 0.05 * sz.c ];
% gaussz1 = make_gaussian( sz.sz, mu, stdev );

for pp = 1 : probe.scpm.N

    %probez0pad( :, :, pp ) = padarray( probez0( :, :, pp ), z0pad );
    probez0pad( :, :, pp ) = zeropadarray( probez0( :, :, pp ), z0pad );
    
end

% for pp = 1 : probe.scpm.N
% 
%     figure;
%     subplot(121); imagescHSV( log10(1 + 1e0 * abs( probez0pad( :, :, pp ) )) .* exp(1 * 1i * angle( probez0pad( :, :, pp ) ))); daspect([1,1,1]);
%     subplot(122); imagesc( log10(1 + 1e0 * abs( probez0pad( :, :, pp ) )) ); daspect([1,1,1]); colormap(expt.cm.blkgrn);
%     
% end
% 
% close all;

%==================================================================================================

% now that we've ensured that we can use spectrum propagation to propagate from the beam 
% defining aperture plane z0 to the sample plane z2, perform this wavefield propagation.

probez2 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

% fov of sample plane z2 is the same as the focal plane z1 since we use spectrum propagation
csys.z2.Lr = csys.z0.Lr;
csys.z2.Lc = csys.z0.Lc;
csys.z2.dr = csys.z0.Lr / sz.r;
csys.z2.dc = csys.z0.Lc / sz.c;

wfp = struct;
wfp.lambda = expt.lambda; 
wfp.zif = +1 * csys.z02; 
wfp.zi = csys.z0;

for pp = 1 : probe.scpm.N
    
    % propagate from to sample plane z2 using spectrum propagation:
    [ probez2( :, :, pp ) ] = fresnelpropspectrum( probez0pad( :, :, pp ), wfp );
    
end

% for pp = 1 : probe.scpm.N
%     
%     figure; 
%     subplot(121); imagescHSV( log10(1 + 10^2 * abs( probez2( :, :, pp ) )) .* exp(1i*angle(probez2( :, :, pp ))) ); daspect([1 1 1])
%     subplot(122); imagesc( log10(1 + 10^2 * abs( probez2( :, :, pp ) ))); daspect([1 1 1]); colormap(expt.cm.blkgrn);
% 
% end
% 
% close all;

% pm = 3;
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^1 * abs(probez2( :, :, pm ))) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(probez2( :, :, pm )) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz02\_sp')        

%==================================================================================================

probe.P = probez2;

for ii = 1 : 1
    
    % orthogonalize the probe modes:
    [ probe.P ] = orthog_modes_eigendecomp( probez2, sz.sz, probe.scpm.N );
    
    % ensure modes have the desired occupancy:
    for pp = 1 : probe.scpm.N
        probe.P( :, :, pp ) = probe.P( :, :, pp ) * expt.Nphotons * probe.scpm.occ( pp ) / norm( probe.P( :, :, pp ), 'fro' );
    end
    
%     [   norm( probe.P( :, :, 1 ), 'fro' ), ...
%         norm( probe.P( :, :, 2 ), 'fro' ), ...
%         norm( probe.P( :, :, 1 ), 'fro' ), ...
%         norm( probe.P( :, :, 1 ), 'fro' ) ]
    
end

% for pp = 1 : probe.scpm.N
%     
%     figure; 
%     subplot(221); imagescHSV( log10(1 + 10^0 * abs( probez2( :, :, pp ) )) .* exp(1i*angle( probez2( :, :, pp ))) ); daspect([1 1 1])
%     subplot(222); imagesc( log10(1 + 10^0 * abs( probez2( :, :, pp ) ))); daspect([1 1 1]); colormap jet
%     subplot(223); imagescHSV( log10(1 + 10^0 * abs( probe.P( :, :, pp ) )) .* exp(1i*angle( probe.P( :, :, pp ))) ); daspect([1 1 1])
%     subplot(224); imagesc( log10(1 + 10^0 * abs( probe.P( :, :, pp ) ))); daspect([1 1 1]); colormap jet
%  
% end
%    
% close all;

return

%==================================================================================================












% 
% 
% % try out z2 to z3 using spectrum propagation:
% 
% z2pad = round( 1.0 * [ 256, 256 ] );
% 
% zcrit2.r = csys.z2.dr ^ 2 * ( 2 * z2pad( 1 )) / expt.lambda;
% zcrit2.c = csys.z2.dc ^ 2 * ( 2 * z2pad( 2 )) / expt.lambda;
% 
% fprintf( [ '\n', num2str( 1e3 * [ zcrit2.r, zcrit2.c ], 'zcrit2 = ( %.3f, %.3f ) mm' ), ...
%                  num2str( 1e3 * csys.z02, ', z02 = %.3f mm'), ...
%                  num2str( 1e3 * csys.z23, ', z23 = %.3f mm'), ...
%                  num2str( 1e3 * csys.z03, ', z03 = %.3f mm'), ' \n\n' ])
%              
%              
%              
% pm = 3;
% 
% % fov of sample plane z2 is the same as the focal plane z1 since we use spectrum propagation
% csys.z3.Lr = csys.z2.Lr;
% csys.z3.Lc = csys.z2.Lc;
% csys.z3.dr = csys.z2.Lr / size( probez2( :, :, pm ), 1 );
% csys.z3.dc = csys.z2.Lc / size( probez2( :, :, pm ), 1 );
% 
% wfp = struct;
% wfp.lambda = expt.lambda; 
% wfp.zif = +1 * csys.z23; 
% wfp.zi = csys.z3;
% 
%   
% % propagate from to sample plane z2 using spectrum propagation:
% [ probez3_z2z3_SP ] = fresnelpropspectrum( probez2( :, :, pm ), wfp );
% 
% 
% 
% pm = 3;
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^1 * abs(probez3_z2z3_SP)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(probez3_z2z3_SP) ); colormap( ax3, hsv ); daspect([1 1 1]); title('probez3\_z2z3\_SP')        
% 
% 
% 
% 
% 
% %================================================

% also try out z0 to z3 using spectrum propagation

probez3 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

% fov of sample plane z2 is the same as the focal plane z1 since we use spectrum propagation
csys.z3.Lr = csys.z0.Lr;
csys.z3.Lc = csys.z0.Lc;
csys.z3.dr = csys.z0.Lr / sz.r;
csys.z3.dc = csys.z0.Lc / sz.c;

wfp = struct;
wfp.lambda = expt.lambda; 
wfp.zif = +1 * csys.z03; 
wfp.zi = csys.z0;

for pp = 1 : probe.scpm.N
    
    % propagate from to sample plane z2 using spectrum propagation:
    [ probez3( :, :, pp ) ] = fresnelpropspectrum( probez0pad( :, :, pp ), wfp );
    
end

% for pp = 1 : probe.scpm.N
%     
%     figure; 
%     subplot(121); imagescHSV( log10(1 + 10^0 * abs( probez3( :, :, pp ) )) .* exp(1i*angle(probez3( :, :, pp ))) ); daspect([1 1 1])
%     subplot(122); imagesc( log10(1 + 10^0 * abs( probez3( :, :, pp ) ))); daspect([1 1 1]); colormap jet
% 
% end
% 
% close all;

pm = 3;
figure; 
ax2 = subplot(121); imagesc( log10(1+ 10^1 * abs(probez3( :, :, pm ))) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
ax3 = subplot(122); imagesc( angle(probez3( :, :, pm )) ); colormap( ax3, hsv ); daspect([1 1 1]); title('probez3\_z0z3\_SP')        



%==================================================================================================

% propagate wavefield z0 to z3 using direct method

pm = 3;

%===============================

% z0 to z3:

phiz0 = probez0pad( :, :, pm );

csys03.z3.Lr = single( expt.lambda * csys.z03 / csys.z0.dr ); 
csys03.z3.Lc = single( expt.lambda * csys.z03 / csys.z0.dc ); 
csys03.z3.dr = csys03.z3.Lr / size( phiz0, 1 );
csys03.z3.dc = csys03.z3.Lc / size( phiz0, 2 );

wfp.lambda = expt.lambda; 
wfp.zif = csys.z03; 
wfp.zi = csys.z0;
wfp.zf = csys03.z3;
wfp.intcurv = 1;
wfp.extcurv = 1;
wfp.dir = 'forward';
    phiz3_z03D = fresnelpropdirect( phiz0, wfp );
    
phiz3_z03D = -1i * phiz3_z03D;
    
slc03.r = round( ( csys03.z3.Lr / csys.z3.Lr ) * size( phiz3_z03D, 1 ));
slc03.r = ceil( ( slc03.r - 2 ) / 2 ) * 2 + 2;
slc03.c = round( ( csys03.z3.Lc / csys.z3.Lc ) * size( phiz3_z03D, 2 ));
slc03.c = ceil( ( slc03.c - 2 ) / 2 ) * 2 + 2;
probez3_tr = truncatearray( probez3( :, :, pm ), [ slc03.r, slc03.c ] );


    
figure; 
ax2 = subplot(121); imagesc( log10(1+ 10^2 * abs(probez3_tr)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
ax3 = subplot(122); imagesc( angle(probez3_tr) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz3\_z03 SP Trunc') 
    
    
figure; 
ax2 = subplot(121); imagesc( log10(1+ 10^2 * abs(phiz3_z03D)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
ax3 = subplot(122); imagesc( angle(phiz3_z03D) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz3\_z03 Direct')



%===============================

% z2 to z3 USING DIRECT METHOD !!!:

phiz2 = probez2( :, :, pm );

%========

% z2pad = 0 * [ 256, 256 ];
% phiz2 = padarray( phiz2, z2pad );

%========

% for the direct method of evaluating the fresnel integral, we need to be able to sample the
% quadratic phase curvature terms adequately. the outside curvature term limits the propagation 
% distance to short distances while the inside term limits it to longer distances.


% rs = 1.7750;
rs = 4;

% zcrit2.r = csys.z2.dr ^ 2 * size( phiz2, 1 ) / expt.lambda;
% zcrit2.c = csys.z2.dc ^ 2 * size( phiz2, 2 ) / expt.lambda;
zcrit2.r = csys.z2.Lr ^ 2 / ( rs * size( phiz2, 1 ) * expt.lambda );
zcrit2.c = csys.z2.Lc ^ 2 / ( rs * size( phiz2, 2 ) * expt.lambda );

fprintf( '\nCheck to see if we can use Spectrum Propagation to get from z0 to z2:');
fprintf( [ '\n', num2str( 1e3 * [ zcrit2.r, zcrit2.c ], 'zcrit2 = ( %.3f, %.3f ) mm' ), ...
                 num2str( 1e3 * csys.z02, ', z02 = %.3f mm'), ...
                 num2str( 1e3 * csys.z23, ', z23 = %.3f mm'), ...
                 num2str( 1e3 * csys.z03, ', z03 = %.3f mm'), ' \n\n' ])
   
%========

% % fft-zeropad interpolation
% rs = 3;
% sz = size( phiz2 );
% tmp0 = fftshift( fft2( fftshift( phiz2 ))) / sqrt( numel( phiz2 ));
% % tmp0 = padarray( tmp0, round( rs * 0.5 * [sz(1), sz(2)] ));
% tmp0 = zeropadarray( tmp0, round( rs * 0.5 * [sz(1), sz(2)] ));
% phiz2 = fftshift( ifft2( fftshift( tmp0 )));

%========

% phiz2 = imresize( phiz2, rs, 'lanczos3' );     % box, triangle, bilinear, nearest
phiz2 = imresize2D( phiz2, rs * size( phiz2 ));

%========

% tmp = fftshift( fft2( fftshift( phiz2pad )));
% tmp = truncatearray( tmp, round( 0.5 * size( tmp )) );
% phiz2pad = fftshift( ifft2( fftshift( tmp )));

%========

% zcrit2.r = csys.z2.dr ^ 2 * ( 2 * z2pad( 1 )) / expt.lambda;
% zcrit2.c = csys.z2.dc ^ 2 * ( 2 * z2pad( 2 )) / expt.lambda;
% 
% fprintf( [ '\n', num2str( 1e3 * [ zcrit2.r, zcrit2.c ], 'zcrit2 = ( %.3f, %.3f ) mm' ), ...
%                  num2str( 1e3 * csys.z02, ', z02 = %.3f mm'), ...
%                  num2str( 1e3 * csys.z03, ', z03 = %.3f mm'),' \n\n' ]);
             
%===============================

csys23.z2.dr = csys.z2.dr * ( 1 / rs );
csys23.z2.dc = csys.z2.dc * ( 1 / rs );

csys23.z3.Lr = single( expt.lambda * csys.z23 / csys23.z2.dr ); 
csys23.z3.Lc = single( expt.lambda * csys.z23 / csys23.z2.dc ); 
csys23.z3.dr = csys23.z3.Lr / size( phiz2, 1 );
csys23.z3.dc = csys23.z3.Lc / size( phiz2, 2 );
           
%===============================

wfp.lambda = expt.lambda; 
wfp.zif = csys.z23; 
wfp.zi = csys23.z2;
wfp.zf = csys23.z3;
wfp.intcurv = 1;
wfp.extcurv = 1;
wfp.dir = 'forward';
    phiz3_z23 = fresnelpropdirect( phiz2, wfp );

phiz3_z23 = -1i * phiz3_z23;
    
%===============================

if csys03.z3.Lr < csys23.z3.Lr
    
    slc23.r = round( ( csys03.z3.Lr / csys23.z3.Lr ) * size( phiz3_z23, 1 ));
    slc23.r = ceil( ( slc23.r - 2 ) / 2 ) * 2 + 2;
    phiz3_z23 = truncatearray( phiz3_z23, [ slc23.r, size( phiz3_z23, 2 ) ] );
    
end

if csys03.z3.Lc < csys23.z3.Lc
    
    slc23.c = round( ( csys03.z3.Lc / csys23.z3.Lc ) * size( phiz3_z23, 2 ));
    slc23.c = ceil( ( slc23.c - 2 ) / 2 ) * 2 + 2;
    phiz3_z23 = truncatearray( phiz3_z23, [ size( phiz3_z23, 1 ), slc23.c ] );
    
end

%===============================

figure; 
ax2 = subplot(121); imagesc( log10(1+ 10^2 * abs(phiz3_z23)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
ax3 = subplot(122); imagesc( angle(phiz3_z23) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz3\_z23')

5;

%===============================

% % now, take phiz3_z23, see if we can down interpolate, truncate, or whatever and get back probez2
% 
% rs = 0.5;
% phiz3_z23_trnc = imresize( phiz3_z23, rs, 'lanczos3' );     % box, triangle, bilinear, nearest
% 
% csys23.z3.dr = csys23.z3.dr * 2;
% csys23.z3.dc = csys23.z3.dc * 2;
% % csys23.z2.Lr = csys23.z2.Lr / 2; 
% % csys23.z2.Lc = csys23.z2.Lc / 2; 
% 
% 
% wfp.lambda = expt.lambda; 
% wfp.zif = -csys.z23; 
% wfp.zi = csys23.z3;
% wfp.zf = csys23.z2;
% wfp.intcurv = 1;
% wfp.extcurv = 1;
% wfp.dir = 'backward';
%     phiz2_z32 = fresnelpropdirect( phiz3_z23_trnc, wfp );
% 
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^1 * abs(phiz2)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(phiz2) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz2')
% 
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^1 * abs(phiz2_z32)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(phiz2_z32) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz2\_z32')

%}

%==================================================================================================

% only use the probe where it's "not too small":

% for pp = 1 : probe.scpm.N
%     
%     tmp0 = abs( probez2( :, :, pp ));
%     tmp1 = max( tmp0( : ));
%     
%     probez2( :, :, pp ) = probez2( :, :, pp ) .* ( tmp0 > 0.01 * tmp1 );
%     
% end
% 
% % for pp = 1 : probe.scpm.N
% %     
% %     figure; 
% %     subplot(131); imagescHSV( log10(1 + 10^0 * abs( probez2( :, :, pp ) )) .* exp(1i*angle(probez2( :, :, pp ))) ); daspect([1 1 1])
% %     ax2 = subplot( 132 ); imagesc( log10(1 + 10^0 * abs( probez2( :, :, pp ) ))); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% %     ax3 = subplot( 133 ); imagesc( abs( probez2( :, :, pp )) > 0.005 * max(max( abs(probez2( :, :, pp ) ))) ); colormap( ax3, gray ); daspect([1 1 1])
% %     
% % end
% % 
% % close all;

%==================================================================================================
