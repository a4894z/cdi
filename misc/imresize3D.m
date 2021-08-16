function [ newimg ] = imresize3D( img, out_prob_sz, iterp_type )

% nearest, linear, cubic, spline, makima
if ~exist( iterp_type,'var' ), iterp_type = 'linear'; end
 
nh = round( out_prob_sz( 1 ));
nw = round( out_prob_sz( 2 ));
nd = round( out_prob_sz( 3 ));

sz = size( img );
ht_scale  = sz( 1 ) / nh;
wid_scale = sz( 2 ) / nw;
dep_scale = sz( 3 ) / nd;

newimg  = interp3( img, ( 1 : 1 : nw )  * wid_scale, ...
                        ( 1 : 1 : nh )  * ht_scale, ...
                        ( 1 : 1 : nd )' * dep_scale, ...
                        iterp_type, 0 );
