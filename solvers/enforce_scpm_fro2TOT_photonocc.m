function [ probe, fro2TOT, scpm_occ ] = enforce_scpm_fro2TOT_photonocc( probe, fro2TOT, scpm_occ )

Nmodes = size( probe, 3 );

fro2 = squeeze( sum( sum( abs( probe ) .^2 )));

%==================

if isempty( fro2TOT )
    
    fro2TOT = sum( fro2 );

end

%==================

is_empty_occ = isempty( scpm_occ );

if Nmodes ~= length( scpm_occ ) && ~is_empty_occ
    
    error('Number of scpm and occupancies do not agree, exiting...');
    
end
    
if is_empty_occ
    
    scpm_occ = fro2 / fro2TOT; 
    
    scpm_occ = sort( scpm_occ / norm( scpm_occ, 1 ));     % make sure the mode occupancy adds up to 1.0
    
end

%==================

% if isempty( fro2TOT ),  fro2TOT  = sum( fro2 ); end

for pp = 1 : Nmodes
    
    probe( :, :, pp ) = probe( :, :, pp ) * sqrt( scpm_occ( pp ) * fro2TOT / fro2( pp ));

end

