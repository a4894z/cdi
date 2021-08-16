function writeimage( p1, opt )

if ~isfield( opt, 'imsz' ),     opt.imsz = [ 512,512 ]; end
if ~isfield( opt, 'phsoff' ),   opt.phsoff = 0; end
if ~isfield( opt, 'logabs' ),   opt.logabs = false; end
if ~isfield( opt, 'bright' ),   opt.bright = 1; end

p1 = p1 / max( abs( p1( : )));

temp1 = complex2hsv( p1, opt.logabs, opt.phsoff, 'no_unwrap' );
temp1( :, :, 3 ) = opt.bright * temp1( :, :, 3 );

temp1 = hsv2rgb( temp1 );

temp2 = imresize( temp1, round( opt.imsz ), 'box' );

%===============================

% 
% sz = size( p1 );
% 
% nh = round( opt.imsz( 1 ));
% nw = round( opt.imsz( 2 ));
% 
% sclR = sz( 1 ) / nh;
% sclC = sz( 2 ) / nw;
% 
% 
% % [ Xi, Yi, Zi ] = meshgrid( (1 : 1 : nw) * sclC, (1 : 1 : nh) * sclR, (1 : 3) );
% %nearest, linear, cubic, or spline
% % temp2  = interp3( temp1, Xi, Yi, Zi, 'nearest', 0 );
% % imwrite( temp2, [opt.name,'.ppm'], 'ppm' );
% 
% %nearest, linear, cubic, or spline
% temp2  = interp3( temp1, (1 : 1 : nw) * sclC, (1 : 1 : nh)' * sclR, (1 : 3), 'nearest', 0 );

%===============================




if strcmp( opt.name( end - 3 : end ), '.ppm' )
  
  imwrite( temp2, opt.name, 'ppm' );

elseif  strcmp( opt.name( end - 3 : end ), '.jpg' )
  
  imwrite( temp2, opt.name, 'jpg', 'Quality', 100 );
  
elseif strcmp( opt.name( end - 3 : end ), '.png' )
  
  ppi = 300;
  dpm = fix( ( ppi / 2.54 ) * 100 );
  imwrite( temp2, opt.name, 'png', 'ResolutionUnit', 'meter', 'XResolution', dpm, 'YResolution', dpm );
  
else
  
  error(' invalid image format you asshat! ')  
  
end

%--------------------------------------------------------------------------------------------------


%{

function wrte_CmplxArr( plot_me, scalefactor, opt.name, wrte_format, wrte_dir, opt.phsoff, opt.logabs, num_iter, quiet )

% wrte_CmplxArr( plot_me, opt.name, opt.logabs, opt.phsoff, wrte_dir, num_iter, quiet )
%
% Takes a complex valued array, and represents it in an HSV colorspace,
% i.e. phase --> hue and magnitude --> value. The saturation component
% is set to 1. The resultant HSV image is then written to a jpf file.
%
% Inputs: 
%
% plot_me   -->   a complex valued 2D array.
%
% opt.name -->   a string which identifies what you want to call the
%                 saved jpg.
%
% wrte_format -->  if 'jpg', write jpg. if 'ppm', write ppm. if 'mat', write mat
%
% wrte_dir -->    a string which identifies the directory you want to 
%                 save the jpg in.
%
% opt.phsoff -->   a real valued scalar, used to give plot_me a global 
%                 phase shift, which in turn will change the hue of the 
%                 image.
%
% opt.logabs --> when we defined the transmission function using some 
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

if ~exist('opt.logabs','var'), opt.logabs = 0; end
if ~exist('opt.phsoff','var'), opt.phsoff = 0; end
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

hsv_img = complex2hsv( plot_me(slcY,slcX), opt.logabs, opt.phsoff );

sz = size(hsv_img);

r_op( 1:sz(1), 1:sz(2), 1:sz(3) ) = hsv2rgb( hsv_img );

%--------------------------------------------------------------------------------------------------

if strcmp(wrte_format,'jpg'),
  
  imwrite(r_op,[wrte_dir, opt.name, num2str(num_iter,'%0.5d'),'.jpg'], 'jpg', 'Quality', 100);
  
elseif strcmp(wrte_format,'ppm'),
  
  imwrite(uint16(65535*r_op),[wrte_dir, opt.name, num2str(num_iter,'%0.5d'),'.ppm'], 'ppm');
  
elseif strcmp(wrte_format,'mat'),
  
  save([wrte_dir, opt.name, num2str(num_iter,'%0.5d')],'plot_me')
  
end
%}



%{
if strcmp( file_opt.name( end - 3 : end ), '.ppm' ) 

  tru.smpl_trans = image2transmission( file_opt.name, tru.abs_min_max, tru.phs_min_max );

elseif strcmp( file_opt.name( end - 3 : end ), '.mat' )

  temp1 = load( file_opt.name );
  tru.smpl_trans = temp1.tru_sampl_trans;

end
%}
%--------------------------------------------------------------------------------------------------
