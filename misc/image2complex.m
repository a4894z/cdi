function [ sample ] = image2complex( file_name, mag, phs )

% 
% From an image, define the sample transmission function using the hue and val 
% components of the HSV representation of the image
% 
% Inputs:
%   file_name  A string with the name of the image to be loaded into Matlab
% 
% Outputs:
%   sample    An array containing the sample transmission function. Dimensions 
%                 are dependent on the image being loaded


%load the image, convert from rgb to hsv
test_img = rgb2hsv( imread( file_name ));

% assign the val component to the magnitude (is scaled between 0 and 1):
test_img_abs = single( test_img( :, :, mag ));

%assign the hue component to the phase (is scaled between 0 and 1):
test_img_phs = single( test_img( :, :, phs ));

sample = test_img_abs .* exp( 1i * test_img_phs );

return



%{

%==================================================================================================

% abs_min_max = [ 0.1, 0.9 ];
% phs_min_max = [ -0.7, -0.2 ] * pi;

% CHECK IF SIZE IS 2 VECTOR
f1 = ~isreal( abs_min_max ) || ~isvector( abs_min_max );
f1 = f1 || abs_min_max(1) < 0 || abs_min_max(2) < 0;
f1 = f1 || (abs_min_max(1) > abs_min_max(2));

if f1, error('improper bounds on abs value of transmission function'); end

f1 = ~isreal( phs_min_max ) || ~isvector( phs_min_max );
f1 = f1 || (abs_min_max(1) > abs_min_max(2));

if f1, error('improper bounds on phase of transmission function'); end


% if ~exist('phs_min_max','var'), phs_min_max = [ -pi + 1e-5, pi - 1e-5 ]; end
% if ~exist('abs_min_max','var'), abs_min_max = [ 0.1, 0.9 ]; end
if ~exist('take_exp_of','var'), take_exp_of = 'no'; end


%==================================================================================================

%load the image, convert from rgb to hsv
test_img = rgb2hsv( imread( file_name ));

%==================================================================================================

% assign the val component to the magnitude:
test_img_abs = single( test_img( :, :, 3 ));

%=======================


% N.Nr = size( test_img_abs, 1 );
% N.Nc = size( test_img_abs, 2 );
% 
% [ test_img_abs, ~ ] = lpf_gauss( test_img_abs, N, N.Nc * 0.1, N.Nr * 0.1 );
% 
% test_img_abs = test_img_abs / max(max( test_img_abs ));


%=======================

% play with modulus contrast:

% test_img_abs = abs( 1 - test_img_abs );
% test_img_abs = test_img_abs / max(max( test_img_abs )) - 0.0;

% test_img_abs = im2bw( test_img_abs );


% test_img_abs = imresize(imresize( test_img_abs, [480,480] ), size(test_img_abs) );
% test_img_abs = test_img_abs / max(max( test_img_abs ));

%increase contrast:
% test_img_abs = exp( 4 * test_img_abs );
% test_img_abs = test_img_abs / max(max( test_img_abs ));


% for ii = 1 : 2
%   test_img_abs = imsharpen( test_img_abs,'Radius',1,'Amount',1.5, 'Threshold', 1.0);
% end

% % 
% abs_min_max = [0.5, 0.4];

%=======================

if min( test_img_abs(:) ) ~= max( test_img_abs(:) )
  
  test_img_abs = test_img_abs - min(min( test_img_abs )); 
  
  test_img_abs = test_img_abs * ( abs_min_max(2) - abs_min_max(1) ) / (1e-7 + max(max(abs(test_img_abs)))) + abs_min_max(1);

end

%==================================================================================================
%phs_min_max = [pi/6, 3*pi/4];

%assign the hue component to the phase:
test_img_phs = single(test_img(:,:,1));
%test_img_phs = test_img_abs;


% N.Nr = size( test_img_phs, 1 );
% N.Nc = size( test_img_phs, 2 );
% 
% [ test_img_phs, ~ ] = lpf_gauss( test_img_phs, N, N.Nc * 0.05, N.Nr * 0.05 );
% 
% test_img_phs = test_img_phs / max(max( test_img_phs ));



if any(test_img_phs(:)) % test if read in image is greyscale (no hue component), if so, skip phase definition (is zero)

  test_img_phs = test_img_phs - min( test_img_phs(:) );
  test_img_phs = test_img_phs / (1e-7 + max(abs(test_img_phs(:))));
    
  test_img_phs = phase_unwrap( angle(exp(1i * 2 * pi * test_img_phs )) );
  
  test_img_phs = test_img_phs - min( test_img_phs(:) );
  test_img_phs = test_img_phs / (1e-7 + max(abs(test_img_phs(:))));
  
  test_img_phs = test_img_phs * ( phs_min_max(2) - phs_min_max(1) ) / (1e-7 + max(abs(test_img_phs(:)))) + phs_min_max(1);
  %test_img_phs = 0.8 * pi * 2 * ( test_img_phs - 0.0 );

end

%==================================================================================================

%define the sample transmission function:

if strcmp( take_exp_of, 'yes' )
  
  sample = single( exp( -1.0 * ( 0.0 + 1 * test_img_abs )) .* exp( 1i * test_img_phs + 1i * 2 * pi * phs_offset ));
  
else
  
  sample = single( ( 0.0 + 1 * test_img_abs ) .* exp( 1i * test_img_phs + 1i * 2 * pi * phs_offset ));

end


%sample = single( exp( -1.0 * (0.0 + 1*test_img_abs)) .* exp( 1i * test_img_phs ) );
% sample = single( (0.0 + 1*test_img_abs) .* exp( 1i * test_img_phs ) );

%==================================================================================================

%clear variables not needed anymore:
clear('test_*');

%}
