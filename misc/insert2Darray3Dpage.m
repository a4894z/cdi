function [ A3d ] = insert2Darray3Dpage( A3d0, B2d, page )



% A = reshape( 1 : ( sr * sc * sp ), sr, sc, sp );
% 
% B = reshape( 2000 + ( 1 : sr * sc )', [ sr, sc ] );
% 
% page = 6;   



if page > ( size( A3d0, 3 ) + 1 ), error('Can only use a insert location of up to 3d array page size + 1'); end
if page < 1, error('insert location must be larger than 1'); end

A3d = cat( 3, A3d0( :, :, 1 : page - 1 ), B2d, A3d0( :, :,  page : end ));
