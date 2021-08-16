function write_complex2hsv_image1x3( p1, p2, p3, wrte_name, out_prob_sz, log_absrptn, globl_phs, brightness_scaling )

%--------------------------------------------------------------------------------------------------

if ~exist('globl_phs','var'), globl_phs = 0; end
if ~exist('log_absrptn','var'), log_absrptn = 0; end
if ~exist('brightness_scaling','var'), brightness_scaling = 1; end

%--------------------------------------------------------------------------------------------------

p1 = p1 / max(abs(p1(:)));
p2 = p2 / max(abs(p2(:)));
p3 = p3 / max(abs(p3(:)));

%--------------------------------------------------------------------------------------------------

szR = size(p1,1);
szC = size(p1,2);

% p2 = im_resize2d( p2, [szR, szC], 'nearest' );
% p3 = im_resize2d( p3, [szR, szC], 'nearest' );
p2 = imresize(p2, [szR, szC]);
p3 = imresize(p3, [szR, szC]);

%--------------------------------------------------------------------------------------------------

p_123( 1 : szR, 1 : szC ) = p1;

p_123( 1 : szR, (szC + 1) : (szC + szC) ) = p2;

p_123( 1 : szR, (szC + 1 + szC) : (szC + szC + szC) ) = p3;

bullseye_color = pi * 1;
bullseye_bright = 0.2;
p_123( szR/2 + 1, 1 : (szC + szC + szC) ) = bullseye_bright * exp(1i * bullseye_color); 
p_123( 1 : szR, szC/2 + 1 ) = bullseye_bright * exp(1i * bullseye_color);
p_123( 1 : szR, szC + szC/2 + 1 ) = bullseye_bright * exp(1i * bullseye_color);
p_123( 1 : szR, 2*szC + szC/2 + 1 ) = bullseye_bright * exp(1i * bullseye_color);

%--------------------------------------------------------------------------------------------------

temp1 = complex2hsv( p_123, log_absrptn, globl_phs, 'no_unwrap' );
temp1(:,:,3) = brightness_scaling * temp1(:,:,3);
temp1 = hsv2rgb( temp1 );

%--------------------------------------------------------------------------------------------------

nh = round( out_prob_sz(1) );
nw = round( 3 * out_prob_sz(2) );

%--------------------------------------------------------------------------------------------------

ht_scale = szR / nh;
wid_scale = 3 * szC / nw;

%--------------------------------------------------------------------------------------------------

% [ Xi, Yi, Zi ] = meshgrid( (1 : 1 : nw) * wid_scale, (1 : 1 : nh) * ht_scale, (1 : 3) );
%nearest, linear, cubic, or spline
% temp2  = interp3( temp1, Xi, Yi, Zi, 'nearest', 0 );
% imwrite( temp2, [wrte_name,'.ppm'], 'ppm' );

%nearest, linear, cubic, or spline
temp2  = interp3( temp1, (1 : 1 : nw) * wid_scale, (1 : 1 : nh)' * ht_scale, (1 : 3), 'nearest', 0 );

%--------------------------------------------------------------------------------------------------

if strcmp( wrte_name( end - 3 : end ), '.ppm' ),
  
  imwrite( temp2, wrte_name, 'ppm' );

elseif  strcmp( wrte_name( end - 3 : end ), '.jpg' ),
  
  imwrite( temp2, wrte_name, 'jpg', 'Quality', 75 );
  
elseif strcmp( wrte_name( end - 3 : end ), '.png' ),
  
  ppi = 600;
  dpm = fix( ( ppi / 2.54 ) * 100 );
  imwrite( temp2, wrte_name, 'png', 'ResolutionUnit', 'meter', 'XResolution', dpm, 'YResolution', dpm );
  
else
  
  error(' invalid image format you asshat! ')
  
end

%--------------------------------------------------------------------------------------------------



%imwrite( temp2, [wrte_name,'.jpg'], 'jpg', 'Quality', 75 );
%imwrite( temp2, [wrte_name,'.ppm'], 'ppm' );



%imwrite( temp2, [wrte_name,'.jpg'], 'jpg', 'Quality', 100 );


% temp2 = imresize( temp1 ,[512 3*512], 'box' );
% imwrite( temp2, wrte_name, 'jpg', 'Quality', 100 );




%{
p_123(1 : szR, 1 : szC, 1:3 ) = complex2hsv( p1, 0, 0 );

p_123( 1 : szR, (szC + 1) : (szC + szC), 1:3 ) = complex2hsv( p2, 0, 0 );

p_123( 1 : szR, (szC + 1 + szC) : (szC + szC + szC), 1:3 ) = complex2hsv( p3, 0, 0 );
%}






end




