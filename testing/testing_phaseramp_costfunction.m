%

%==================================================================================================
% create an exit wave:

N   = [ 128, 128 ];
len = 0.3 * N;

[ exwv ] = make_rectangle( N, len );

%==================================================================================================
% create a phase ramp in the exit wave:

% phase ramp slope
prs = [ 1.0, 3.0 ];

[ c, r ] = meshgrid( ( -0.5 * N( 2 ) + 1 ) : ( 0.5 * N( 2 )),( -0.5 * N( 1 ) + 1 ) : ( 0.5 * N( 1 )) );
c = single( c ); 
r = single( r ); 

V_c = exp( 1i * 2 * pi * prs( 2 ) * c / N( 2 ));
V_r = exp( 1i * 2 * pi * prs( 1 ) * r / N( 1 )); 

V = V_r .* V_c;

exwv = exwv .* V;

F_exwv = fftshift( fft2( fftshift( exwv )) / sqrt( numel( exwv )));

figure; imagescHSV( exwv )

% figure; imagescHSV( log10( 1 + abs( F_exwv )) .* exp( 1i * angle(F_exwv)) );
figure; imagesc( log10( 1 + abs( F_exwv ))  ); colormap jet

figure; imagesc( angle( exwv ))
figure; imagesc( real( exwv ))
figure; imagesc( imag( exwv ))

figure; imagesc( diff(angle( exwv ), 1, 1 ))
figure; imagesc( diff(angle( exwv ), 1, 2 ))

%==================================================================================================










h5info /home/ash/Desktop/Star_fine_mid_out.h5
h5disp('/home/ash/Desktop/Star_fine_mid_out.h5')
h5disp('/home/ash/Desktop/Star_fine_mid_out.h5','/img')



% STRETCHING A MEASUREMENT
FsampleTF = fftshift( fft2( fftshift( expt.sample.TF ))) / expt.sample.sz.sqrt_rc;
% FsampleTF = padarray( FsampleTF, [ 0, 1024 ] );
FsampleTFstretch = imresize( FsampleTF, round( [ expt.sample.sz.r, 1024 ] ), 'box' );
sampleTFstretch = fftshift( ifft2( fftshift( FsampleTFstretch ))) * expt.sample.sz.sqrt_rc;

figure; 
subplot(121); imagesc( log10(1 + abs(FsampleTF))); daspect([1 1 1])
subplot(122); imagesc( log10(1 + abs(FsampleTFstretch))); daspect([1 1 1])

figure; 
subplot(141); imagesc(abs(sampleTFstretch)); daspect([1 1 1])
subplot(142); imagesc(angle(sampleTFstretch)); daspect([1 1 1])
subplot(143); imagescHSV(sampleTFstretch); daspect([1 1 1])
subplot(144); imagescHSV(expt.sample.TF); daspect([1 1 1])