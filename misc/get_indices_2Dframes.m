function [ ind ] = get_indices_2Dframes( spos_rs, sz, vs_r, vs_c )


% ( r, c ) indices for each scan position:
% vsr_s = uint32( transpose( round( vs.r - 1.0 * spos_rs( :, 1 ))));
% vsc_s = uint32( transpose( round( vs.c - 1.0 * spos_rs( :, 2 ))));
% vsr_s = single( transpose( round( vs.r - 1.0 * spos_rs( :, 1 ))));
% vsc_s = single( transpose( round( vs.c - 1.0 * spos_rs( :, 2 ))));
vsr_s = transpose( round( vs_r - 1.0 * spos_rs( :, 1 )));
vsc_s = transpose( round( vs_c - 1.0 * spos_rs( :, 2 )));

% create total, "matching" ( r, c ) indices for each scan position:
vsr = repmat( vsr_s, size( vsc_s, 1 ), 1 );
vsc = repelem( vsc_s, size( vsr_s, 1 ), 1 );

% ( r, c ) indices to linear indices:
% ind = sub2ind( uint32( sz ), vsr, vsc );
ind = sub2ind( sz , vsr, vsc );
% ind = transpose( ind );

% ind = single( ind );
ind = uint32( ind );


% 
% A = [ 5, 4, 2, 6; ...
%       1, 5, 6, 2 ]; % A is a 2 by n matrix, here n is 4
%   
% B = [ 12, 15, 2, 9, 8, 6 ]'; 
% 
% C1 = [B(5), B(4), B(2), B(6); B(1), B(5), B(6), B(2)]; 
% C2 = B(A)

