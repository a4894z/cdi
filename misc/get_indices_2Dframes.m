function [ ind ] = get_indices_2Dframes( spos_rs, sz, vs_r, vs_c )

% ( r, c ) indices for each scan position:
vsr_s = transpose( round( vs_r - 1.0 * spos_rs( :, 1 )));
vsc_s = transpose( round( vs_c - 1.0 * spos_rs( :, 2 )));

% create total, "matching" ( r, c ) indices for each scan position:
vsr = repmat(  vsr_s, size( vsc_s, 1 ), 1 );
vsc = repelem( vsc_s, size( vsr_s, 1 ), 1 );

% ( r, c ) indices to linear indices:
ind = uint32( sub2ind( sz, vsr, vsc ));

