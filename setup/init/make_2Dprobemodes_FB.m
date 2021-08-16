function [ probe, csys, sz ] = make_2Dprobemodes_FB( expt )

% make_2Dprobemodes_F
% make_2Dprobemodes_UF

%==================================================================================================
%
% We need to define a 2D focusing optic transmission function at z0 plane, then propagate that 
% wavefield to the focal plane location at z1. We don't really care about exact dimensions at z0
% for the focusing optic diameter or pixel size, just the rough shape so that we can numerically
% propagate the wavefield at the focus z1 to the sample plane z2.
%
%
%   |           |          |           |
%   |           |          |           |
%   |           |          |           |
%   |           |          |           |
%   |           |          |           |
%   |           |          |           |
%   |           |          |           |
%   z0          z1         z2          z3
%   focusing    focal      sample      measurement
%   optic       plane      plane       plane
%   pupil 
%   function 
%   plane 
%
%==================================================================================================

csys.z01 = 1.00e-3;                  % focusing optic plane z0 to focal plane z1, i.e. the focal length ( in meters )
csys.z12 = 0.01e-3;                  % focal plane z1 at sample plane z2 ( in meters )
csys.z02 = csys.z01 + csys.z12;      % focusing optic plane z0 to sample plane z2 ( in meters )

csys.z23 = 3.00e-0;                  % sample to detector ( in meters )
csys.z13 = csys.z12 + csys.z23;      % focal plane z1 at measurement plane z3 ( in meters )
csys.z03 = csys.z02 + csys.z23;      % focusing optic plane z0 to measurement plane z3 ( in meters )

%================================================

paths.probe = cell( 0 );

% paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'PETN2.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'pinholeA.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'pinholeC.ppm' ];
% % paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'zp6_512px.ppm' ];
% % paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'spiral_zp.ppm' ];
% % paths.probe{ end + 1 } = [ expt.paths.rdata, '/probe/', 'zp6b_512px.ppm' ];

% paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'zp6b_512px.ppm' ];
% paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'PETN2.ppm' ];
paths.probe{ end + 1 } = [ expt.paths.rdata, 'probe/', 'pinholeEc.ppm' ];

%================================================

% (s)patially (c)oherent (p)robe (m)ode occupancies:
probe.scpm.occ = [ 1.0 ];
% probe.scpm.occ = [ 0.025, 0.075, 0.9 ];
% scpm.occ = [ 0.1, 0.15, 0.2, 0.25, 0.3 ];
% scpm.occ = [ 0.02, 0.03, 0.1, 0.15, 0.7 ];
% scpm.occ = random('unif', 0, 1, length( paths.probe ));

% make sure the mode occupancy adds up to 1.0:
probe.scpm.occ = sort( probe.scpm.occ / norm( probe.scpm.occ, 1 ));

% get total number of scpm used to define simulated x-ray beam:
probe.scpm.N = length( probe.scpm.occ );

%================================================

% problem size for probe, exit wave, measurements.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% this might be different than the final problem 
% size since we maybe will need to zero pad in order 
% to allow use of spectrum propagation from any focal
% plane z1 to sample plane z2 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sz.r = round( 768 ); 
sz.c = round( 768 );
sz.sz = [ sz.r, sz.c ]; 
sz.rc = sz.r * sz.c;
sz.sqrt_rc = sqrt( sz.rc );

%==================================================================================================

% pixel size of the focusing optic pupil function plane z0:
csys.z0.dr = 50.0e-9;
csys.z0.dc = 50.0e-9;
      
% focusing optic pupil function dimensions:
PFz0.Lr = 10e-6;  % vertical diameter/side length ( in meters )
PFz0.Lc = 10e-6;  % horizontal diameter/side length ( in meters )

PFz0.r = round( PFz0.Lr / csys.z0.dr );  % rectangular vertical length ( in pixels )
PFz0.c = round( PFz0.Lc / csys.z0.dc );  % rectangular horizontal length ( in pixels )

if PFz0.r > sz.r, error('Must increase the number of ROW pixels in the focusing optic pupil function plane z0.'); end
if PFz0.c > sz.c, error('Must increase the number of COLUMNS pixels in the focusing optic pupil function plane z0.'); end

clear( 'PFz0' );

%==================================================================================================

% we'll use the single-fft2 "direct" method of wavefield propagation 
% to go from the focusing optic plane z0 to the focal plane z1. 
% so, the field of view at the focal plane z1 is:

% fov of focal plane z1 using focusing optic plane z0 parameters.
csys.z1.Lr = single( expt.lambda * csys.z01 / csys.z0.dr ); 
csys.z1.Lc = single( expt.lambda * csys.z01 / csys.z0.dc ); 

% pixel sizes at the focal plane z1.
csys.z1.dr = csys.z1.Lr / sz.r;
csys.z1.dc = csys.z1.Lc / sz.c;

% % fov of focal plane z1 using focusing optic plane z0 parameters.
% csys.z3.Lr = single( expt.lambda * csys.z13 / csys.z1.dr ); 
% csys.z3.Lc = single( expt.lambda * csys.z13 / csys.z1.dc ); 
% fprintf( num2str( 1e3 * [ csys.z3.Lr, csys.z3.Lc ], 'FOV at z3: %.3f mm, %.3f mm' ) );

%==================================================================================================

% create non-orthogonal modes at the focusing optic z0

probez0 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

abs_minmax_z0 = repmat( [ 0, 1 ], probe.scpm.N, 1);
phs_minmax_z0 = repmat( 0.3 * [ -pi, pi ], probe.scpm.N, 1);

for pp = 1 : probe.scpm.N
   
    tmpz0 = image2complex( paths.probe{ pp } );
%     tmpz0 = padarray( tmpz0, [ 100, 100 ] );
    tmpz0 = imresize( tmpz0, sz.sz, 'nearest' );  % box, triangle, bilinear, nearest
    tmpz0 = lpf_gauss( tmpz0, 0.30 * sz.sz ); 

    tmpz0 = modulus_limits_scale( tmpz0, abs_minmax_z0( pp, : ));
    tmpz0 = phase_limits_scale( tmpz0, phs_minmax_z0( pp, : ));
    
%     [ min( abs(tmpz0(:))), max( abs(tmpz0(:))) ]
%     [ min( angle(tmpz0(:))), max( angle(tmpz0(:))) ]

    probez0( :, :, pp ) = tmpz0;
    
end

% for pp = 1 : probe.scpm.N
%    
%     figure; 
%     subplot(121); imagesc( log10(1 + 1e0 * abs( probez0( :, :, pp ) ))); daspect([1 1 1])
%     subplot(122); imagescHSV( log10(1 + 1e0 * abs( probez0( :, :, pp ))) .* exp(1i * angle(probez0( :, :, pp )))); daspect([1 1 1])
% end
% 
% close all;

%==================================================================================================
% numerically propagate the wavefield at the pupil function plane z0 to the focal plane z1

probez1 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

wfp.lambda = expt.lambda; 
wfp.zif = csys.z01; 
wfp.zi = csys.z0;
wfp.zf = csys.z1;
wfp.intcurv = false;
wfp.extcurv = true;
wfp.dir = 'forward';

for pp = 1 : probe.scpm.N
   
    % propagate from pupil function plane z0 to focal plane z1:
    [ probez1( :, :, pp ) ] = fresnelpropdirect( probez0( :, :, pp ), wfp );

end

% for pp = 1 : probe.scpm.N
%     
%     figure;
%     subplot(121); imagesc( log10(1 + 1e0 * abs( probez1( :, :, pp ) )) ); daspect([1 1 1]); colormap jet
%     subplot(122); imagescHSV( log10(1 + 1e0 * abs( probez1( :, :, pp ) )) .* exp(1 * 1i * angle( probez1( :, :, pp ) ))); daspect([1 1 1]);
% 
% end
% 
% close all;
   
%==================================================================================================
% zero pad the wavefield at the focal plane z1, so we could use spectrum propagation to 
% numerically propagate the wavefield at the focal plane z1 to the sample plane z2

% z1pad = round( 0 * [ 256, 256 ] + 0.5 * 56 );
z1pad = round( 0 * [ 256, 256 ] + 0.0 * 24 );

sz.r = sz.r + ( 2 * z1pad( 1 ));
sz.c = sz.c + ( 2 * z1pad( 2 ));
sz.sz = [ sz.r, sz.c ]; 
sz.rc = sz.r * sz.c;
sz.sqrt_rc = sqrt( sz.rc );

csys.z1.Lr = csys.z1.dr * sz.r;
csys.z1.Lc = csys.z1.dc * sz.c;

%===============================

zcrit1.r = csys.z1.dr ^ 2 * ( 2 * z1pad( 1 )) / expt.lambda;
zcrit1.c = csys.z1.dc ^ 2 * ( 2 * z1pad( 2 )) / expt.lambda;

fprintf( [ '\n', num2str( 1e3 * [ zcrit1.r, zcrit1.c ], 'zcrit1 = ( %.3f, %.3f ) mm' ), ...
                 num2str( 1e3 * csys.z12, ', z12 = %.3f mm'), ...
                 num2str( 1e3 * csys.z23, ', z23 = %.3f mm'),' \n' ])

%===============================

probez1pad = zeros( sz.r, sz.c , 'single' );

% mu = [ sz.r / 2 + 1, sz.c / 2 + 1 ];
% stdev = [  0.05 * sz.r, 0.05 * sz.c ];
% gaussz1 = make_gaussian( sz.sz, mu, stdev );

for pp = 1 : probe.scpm.N

    probez1pad( :, :, pp ) = padarray( probez1( :, :, pp ), z1pad );
    
end

% for pp = 1 : probe.scpm.N
%     
%     figure;
%     subplot(121); imagescHSV( log10(1 + 1e0 * abs( probez1pad( :, :, pp ) )) .* exp(1 * 1i * angle( probez1pad( :, :, pp ) ))); daspect([1,1,1]);
%     subplot(122); imagesc( log10(1 + 1e0 * abs( probez1pad( :, :, pp ) )) ); daspect([1,1,1]); colormap jet;
% 
% end
% 
% close all;

%==================================================================================================

% now that we've ensured that we can use spectrum propagation to propagate from 
% the focal plane z1 to the sample plane z2, perform this wavefield propagation.

probez2 = zeros( sz.r, sz.c, probe.scpm.N, 'single' );

% fov of sample plane z2 is the same as the focal plane z1 since we use spectrum propagation
csys.z2.Lr = csys.z1.Lr;
csys.z2.Lc = csys.z1.Lc;
csys.z2.dr = csys.z2.Lr / sz.r;
csys.z2.dc = csys.z2.Lc / sz.c;

wfp = struct;
wfp.lambda = expt.lambda; 
wfp.zif = +1 * csys.z12; 
wfp.zi = csys.z1;

for pp = 1 : probe.scpm.N
    
    % propagate to sample plane using spectrum propagation:
    [ probez2( :, :, pp ) ] = fresnelpropspectrum( probez1pad( :, :, pp ), wfp );

end

for pp = 1 : probe.scpm.N
    
    figure; 
    subplot(131); imagesc( abs( probez2( :, :, pp ) )); daspect([1 1 1]); colormap jet
    subplot(132); imagesc( log10(1 + 10^1 * abs( probez2( :, :, pp ) ))); daspect([1 1 1]); colormap jet
    subplot(133); imagescHSV( log10(1 + 10^1 * abs( probez2( :, :, pp ) )) .* exp(1i*angle(probez2( :, :, pp ))) ); daspect([1 1 1])

    
end

close all;

%==================================================================================================


%==================================================================================================






% 
% 
% % propagate wavefield at z2 to z3 and compare to z0 to z3:
% 
% pm = 3;
% 
% %===============================
% 
% % z1 to z3:
% 
% phiz1 = probez1pad( :, :, pm );
% 
% csys13.z3.Lr = single( expt.lambda * csys.z13 / csys.z1.dr ); 
% csys13.z3.Lc = single( expt.lambda * csys.z13 / csys.z1.dc ); 
% csys13.z3.dr = csys13.z3.Lr / size( phiz1, 1 );
% csys13.z3.dc = csys13.z3.Lc / size( phiz1, 2 );
% 
% %===============================
% 
% % z2 to z3:
% 
% phiz2 = probez2( :, :, pm );
% 
% %========
% 
% z2pad = 0 * [ 256, 256 ];
% phiz2 = padarray( phiz2, z2pad );
% 
% %========
% 
% % % fft-zeropad interpolation
% % sz = size( phiz2 );
% % phiz2 = fftshift( ifft2( fftshift( padarray( fftshift( fft2( fftshift( phiz2 ))), round( 0.5 * [sz(1), sz(2)]) ))));
% 
% %========
% 
% rs = 3;
% 
% phiz2 = imresize( phiz2, rs );     % box, triangle, bilinear, nearest
% 
% %===============================
% 
% 
% % tmp = fftshift( fft2( fftshift( phiz2pad )));
% % tmp = truncatearray( tmp, round( 0.5 * size( tmp )) );
% % phiz2pad = fftshift( ifft2( fftshift( tmp )));
% 
% zcrit2.r = csys.z2.dr ^ 2 * ( 2 * z2pad( 1 )) / expt.lambda;
% zcrit2.c = csys.z2.dc ^ 2 * ( 2 * z2pad( 2 )) / expt.lambda;
% 
% fprintf( [ '\n', num2str( 1e3 * [ zcrit2.r, zcrit2.c ], 'zcrit2 = ( %.3f, %.3f ) mm' ), ...
%                  num2str( 1e3 * csys.z02, ', z02 = %.3f mm'), ...
%                  num2str( 1e3 * csys.z03, ', z03 = %.3f mm'),' \n\n' ]);
% 
% %===============================
%              
% csys23.z2.dr = csys.z2.dr * ( 1 / rs );
% csys23.z2.dc = csys.z2.dc * ( 1 / rs );
% % csys23.z2.Lr = csys23.z2.dr * size( phiz2, 1 );
% % csys23.z2.Lc = csys23.z2.dc * size( phiz2, 2 );
% 
% csys23.z3.Lr = single( expt.lambda * csys.z23 / csys23.z2.dr ); 
% csys23.z3.Lc = single( expt.lambda * csys.z23 / csys23.z2.dc ); 
% csys23.z3.dr = csys23.z3.Lr / size( phiz2, 1 );
% csys23.z3.dc = csys23.z3.Lc / size( phiz2, 2 );
%            
% %===============================
% 
% wfp.lambda = expt.lambda; 
% wfp.zif = csys.z13; 
% wfp.zi = csys.z1;
% wfp.zf = csys13.z3;
% wfp.intcurv = 1;
% wfp.extcurv = 0;
% wfp.dir = 'forward';
%     phiz3_z13 = fresnelpropdirect( phiz1, wfp );
% 
% wfp.lambda = expt.lambda; 
% wfp.zif = csys.z23; 
% wfp.zi = csys23.z2;
% wfp.zf = csys23.z3;
% wfp.intcurv = 1;
% wfp.extcurv = 0;
% wfp.dir = 'forward';
%     phiz3_z23 = fresnelpropdirect( phiz2, wfp );
% 
% %===============================
% 
% % % use if csys13.z3.Lr,c > csys23.z3.Lr,c 
% % slc13.r = round( ( csys23.z3.Lr / csys13.z3.Lr ) * size( phiz3_z13, 1 ));
% % slc13.r = ceil( ( slc13.r - 2 ) / 2 ) * 2 + 2;
% % slc13.c = round( ( csys23.z3.Lc / csys13.z3.Lc ) * size( phiz3_z13, 2 ));
% % slc13.c = ceil( ( slc13.c - 2 ) / 2 ) * 2 + 2;
% % phiz3_z13 = truncatearray( phiz3_z13, [ slc13.r, slc13.c ] );
% 
% % % use if csys13.z3.Lr,c < csys23.z3.Lr,c 
% % slc23.r = round( ( csys13.z3.Lr / csys23.z3.Lr ) * size( phiz3_z23, 1 ));
% % slc23.r = ceil( ( slc23.r - 2 ) / 2 ) * 2 + 2;
% % slc23.c = round( ( csys13.z3.Lc / csys23.z3.Lc ) * size( phiz3_z23, 2 ));
% % slc23.c = ceil( ( slc23.c - 2 ) / 2 ) * 2 + 2;
% % phiz3_z23 = truncatearray( phiz3_z23, [ slc23.r, slc23.c ] );
% 
% %===============================
% 
% figure; imagescHSV( log10(1+ 10^1 * abs(phiz3_z13)) .* exp( 1 * 1i*angle(phiz3_z13)) ); title('phiz3\_z03')
% figure; imagescHSV( log10(1+ 10^1 * abs(phiz3_z23)) .* exp( 1 * 1i*angle(phiz3_z23)) ); title('phiz3\_z23')
% 
% % phiz3_z13 = phiz3_z13 * 1e3 / norm(phiz3_z13, 'fro');
% % phiz3_z23 = phiz3_z23 * 1e3 / norm(phiz3_z23, 'fro');
% 
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^0 * abs(phiz3_z13)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(phiz3_z13) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz3\_z03')
% 
% figure; 
% ax2 = subplot(121); imagesc( log10(1+ 10^0 * abs(phiz3_z23)) ); colormap( ax2, expt.cm.blkgrn ); daspect([1 1 1])
% ax3 = subplot(122); imagesc( angle(phiz3_z23) ); colormap( ax3, hsv ); daspect([1 1 1]); title('phiz3\_z23')
% 
% 
% 
% 
% 
















%==================================================================================================

for ii = 1 : 1
    
    % ensure modes have the desired occupancy:
    for pp = 1 : probe.scpm.N
        probez2( :, :, pp ) = probez2( :, :, pp ) * expt.Nphotons * probe.scpm.occ( pp ) / norm( probez2( :, :, pp ), 'fro' );
    end
    
    % orthogonalize the probe modes:
    [ probe.P ] = orthog_modes_eigendecomp( probez2, size( probez2 ), probe.scpm.N );

    % ensure modes have the desired occupancy:
    for pp = 1 : probe.scpm.N
        probe.P( :, :, pp ) = probe.P( :, :, pp ) * expt.Nphotons * probe.scpm.occ( pp ) / norm( probe.P( :, :, pp ), 'fro' );
    end

end

for pp = 1 : probe.scpm.N
    
    figure; 
    subplot(131); imagesc( abs( probe.P( :, :, pp ))); daspect([1 1 1]); colormap jet
    subplot(132); imagescHSV( log10(1 + 10^0 * abs( probe.P( :, :, pp ) )) .* exp(1i*angle(probe.P( :, :, pp ))) ); daspect([1 1 1])
    subplot(133); imagesc( log10(1 + 10^0 * abs( probe.P( :, :, pp ) ))); daspect([1 1 1]); colormap jet
 
end
   
close all;

%==================================================================================================

% expt = orderfields( expt );
% 
% fov = expt.fov;
% sz = sz;



%==================================================================================================

%{

for ii = 1 : scpm.N

%     figure; 
%     imagesc( abs( probe.P( :, :, ii ))); 
%     daspect([1 1 1]); 
%     colorbar; 
%     colormap jet;
    
    figure; 
    imagescHSV( probe.P( :, :, ii )); 

end

for ii = 1 : scpm.N
    
    m1 = probe.P( :, :, ii ); 
    
    for jj = 1 : scpm.N
        
        m2 = probe.P( :, :, jj );
        
        probe_ip( ii, jj ) = m1(:)' * m2(:);
        
    end

end

figure; 
subplot(121); imagesc( abs(probe_ip ))
subplot(122); imagesc( log10(1+abs( probe_ip )))
5;
close all;

%}




%==================================================================================================

% a = expt.z02 * expt.lambda / dL0.r ^ 2;
% round( (a - 1 ) / 2 ) * 2 + 1            % round to odd
% round( (a - 2 ) / 2 ) * 2 + 2            % round to even

% szpad.r = sz.r;
% szpad.c = sz.c;
% 
% if zcrit.r < csys.z12
%     
%     szpad.r = csys.z12 * expt.lambda / csys.z1.dr ^ 2;
%     szpad.r = ceil( ( szpad.r - 2 ) / 2 ) * 2 + 2;   % round to next larger even number of pixels
%     
%     szpad.r = szpad.r + 256;
% %     szpad.r = 512;
%     
%     fprintf( [ '\nProbe Modes, Unfocused Probe: Need to zeropad ROWS by ', num2str( szpad.r - sz.r, '%d' ), ' pixels !\n' ] );
%     
% end
% 
% if zcrit.c < csys.z12
%     
%     szpad.c = csys.z12 * expt.lambda / csys.z1.dc ^ 2;
%     szpad.c = ceil( ( szpad.c - 2 ) / 2 ) * 2 + 2;   % round to next larger even number of pixels
%     
%     szpad.c = szpad.c + 256;
% %     szpad.c = 512;
%     
%     fprintf( [ '\nProbe Modes, Unfocused Probe: Need to zeropad COLS by ', num2str( szpad.c - sz.c, '%d' ), ' pixels !\n' ] );
% 
% end


