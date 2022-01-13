function [ phi, TFview ] = enforce_2DTPAsposview( P, TF, vs_r, vs_c, spos, shifttype )

TFview = getview_2DsampleTF( TF, vs_r, vs_c, spos, shifttype );

phi = P .* TFview;


% getview_2DsampleTF( TF, vs_r, vs_c, spos, shifttype )




% Nmodes = size( P, 3 );
% phi = P .* repmat( TFview, 1, 1, szP( 3 ));




