function newimg = imresize2D( img, prob_size, iterp_type )

%IM_RESIZE Resize an image using nearest, linear, spline, or bicubic interpolation
%
%          NEWIMG = IM_RESIZE( IMG, PROB_SIZE, INTERP_TYPE ) Given input image IMG,
%          returns a new image NEWIMG of size PROB_SIZE = [ NH, NW ].

%if nargin ~= 3, error('usage: im_resize( image, prob_size, iterp_type )'); end

if ~exist( 'iterp_type', 'var' )
    iterp_type = 'linear';
end
    
nh = prob_size( 1 );
nw = prob_size( 2 );

ht_scale = size( img, 1 ) / nh;
wid_scale = size( img, 2 ) / nw;


% newimg_re = interp2( real(img), (1:nw) * wid_scale, (1:nh)' * ht_scale, 'cubic' );
% newimg_im = interp2( imag(img), (1:nw) * wid_scale, (1:nh)' * ht_scale, 'cubic' );
% newimg = newimg_re + 1i*newimg_im;


newimg = interp2( img, (1:nw) * wid_scale, (1:nh)' * ht_scale, iterp_type, 0 );  % 'linear', 'nearest', 'cubic', 'spline', 'makima'

