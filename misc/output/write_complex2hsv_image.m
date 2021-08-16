function write_complex2hsv_image( p1, wrte_name, out_prob_sz, log_absrptn, globl_phs, brightness_scaling )

%--------------------------------------------------------------------------------------------------

if ~exist('globl_phs','var'), globl_phs = 0; end
if ~exist('log_absrptn','var'), log_absrptn = 0; end
if ~exist('brightness_scaling','var'), brightness_scaling = 1; end

%--------------------------------------------------------------------------------------------------

p1 = p1 / max(abs(p1(:)));

%--------------------------------------------------------------------------------------------------

szR = size(p1,1);
szC = size(p1,2);

%--------------------------------------------------------------------------------------------------

temp1 = complex2hsv( p1, log_absrptn, globl_phs, 'no_unwrap' );

temp1(:,:,3) = brightness_scaling * temp1(:,:,3);

temp1 = hsv2rgb( temp1 );

%--------------------------------------------------------------------------------------------------

nh = round( out_prob_sz(1) );
nw = round( out_prob_sz(2) );

%--------------------------------------------------------------------------------------------------

ht_scale = szR / nh;
wid_scale = szC / nw;

%--------------------------------------------------------------------------------------------------

% [ Xi, Yi, Zi ] = meshgrid( (1 : 1 : nw) * wid_scale, (1 : 1 : nh) * ht_scale, (1 : 3) );
%nearest, linear, cubic, or spline
% temp2  = interp3( temp1, Xi, Yi, Zi, 'nearest', 0 );
% imwrite( temp2, [wrte_name,'.ppm'], 'ppm' );

%nearest, linear, cubic, or spline
temp2  = interp3( temp1, (1 : 1 : nw) * wid_scale, (1 : 1 : nh)' * ht_scale, (1 : 3), 'nearest', 0 );

%--------------------------------------------------------------------------------------------------

if strcmp( wrte_name( end - 3 : end ), '.ppm' )
  
  imwrite( temp2, wrte_name, 'ppm' );

elseif  strcmp( wrte_name( end - 3 : end ), '.jpg' )
  
  imwrite( temp2, wrte_name, 'jpg', 'Quality', 100 );
  
elseif strcmp( wrte_name( end - 3 : end ), '.png' )
  
  ppi = 300;
  dpm = fix( ( ppi / 2.54 ) * 100 );
  imwrite( temp2, wrte_name, 'png', 'ResolutionUnit', 'meter', 'XResolution', dpm, 'YResolution', dpm );
  
else
  
  error(' invalid image format you asshat! ')  
  
end

%--------------------------------------------------------------------------------------------------


%{

function wrte_CmplxArr( plot_me, scalefactor, wrte_name, wrte_format, wrte_dir, globl_phs, log_absrptn, num_iter, quiet )

% wrte_CmplxArr( plot_me, wrte_name, log_absrptn, globl_phs, wrte_dir, num_iter, quiet )
%
% Takes a complex valued array, and represents it in an HSV colorspace,
% i.e. phase --> hue and magnitude --> value. The saturation component
% is set to 1. The resultant HSV image is then written to a jpf file.
%
% Inputs: 
%
% plot_me   -->   a complex valued 2D array.
%
% wrte_name -->   a string which identifies what you want to call the
%                 saved jpg.
%
% wrte_format -->  if 'jpg', write jpg. if 'ppm', write ppm. if 'mat', write mat
%
% wrte_dir -->    a string which identifies the directory you want to 
%                 save the jpg in.
%
% globl_phs -->   a real valued scalar, used to give plot_me a global 
%                 phase shift, which in turn will change the hue of the 
%                 image.
%
% log_absrptn --> when we defined the transmission function using some 
%                 image, we either used the value component of the image
%                 or exp(-1 * value component of the image). 
%                 The exp(-1 ... ) is also what we would realistically
%                 be recovering in an experiment. So, if you set this
%                 value to a real valued scalar of 1, you'll take the 
%                 log of plot_me before setting it to the value component.
%
% num_iter -->    an additional real valued scalar used to identify the         
%                 image.
%
% quiet -->       set to 1 to get text feedback when this function is 
%                 called.
%                 
%
%
% Outputs:
%
% none
%
%

%--------------------------------------------------------------------------------------------------

if ~exist('log_absrptn','var'), log_absrptn = 0; end
if ~exist('globl_phs','var'), globl_phs = 0; end
if ~exist('wrte_dir','var'), wrte_dir = './'; end
if ~exist('num_iter','var'), num_iter = []; end
if ~exist('quiet','var'), quiet = 1; end

%--------------------------------------------------------------------------------------------------

if quiet == 0, fprintf('\n plotting jpg of reconstruction! \n\n'); end

%--------------------------------------------------------------------------------------------------

%rescale what we want to plot:
%scalefactor = 4;
%scalefactor = [512 512];
%plot_me = imresize( plot_me, scalefactor );
plot_me = im_resize2d( plot_me, size(plot_me) *  scalefactor, 'nearest' );

sz = size(plot_me);

%--------------------------------------------------------------------------------------------------

%subset of object array to plot (for if we want to "zoom in" on particular ROI):
slcY = round( (sz(1)/2 - sz(1)/2 + 1) : (sz(1)/2 + sz(1)/2) );
slcX = round( (sz(2)/2 - sz(2)/2 + 1) : (sz(2)/2 + sz(2)/2) );

%--------------------------------------------------------------------------------------------------

hsv_img = complex2hsv( plot_me(slcY,slcX), log_absrptn, globl_phs );

sz = size(hsv_img);

r_op( 1:sz(1), 1:sz(2), 1:sz(3) ) = hsv2rgb( hsv_img );

%--------------------------------------------------------------------------------------------------

if strcmp(wrte_format,'jpg'),
  
  imwrite(r_op,[wrte_dir, wrte_name, num2str(num_iter,'%0.5d'),'.jpg'], 'jpg', 'Quality', 100);
  
elseif strcmp(wrte_format,'ppm'),
  
  imwrite(uint16(65535*r_op),[wrte_dir, wrte_name, num2str(num_iter,'%0.5d'),'.ppm'], 'ppm');
  
elseif strcmp(wrte_format,'mat'),
  
  save([wrte_dir, wrte_name, num2str(num_iter,'%0.5d')],'plot_me')
  
end
%}



%{
if strcmp( file_name( end - 3 : end ), '.ppm' ) 

  tru.smpl_trans = image2transmission( file_name, tru.abs_min_max, tru.phs_min_max );

elseif strcmp( file_name( end - 3 : end ), '.mat' )

  temp1 = load( file_name );
  tru.smpl_trans = temp1.tru_sampl_trans;

end
%}
%--------------------------------------------------------------------------------------------------