function [ Zorth ] = orthog_modes_eigendecomp( Z )

Nmodes = size( Z, 3 );

if Nmodes == 1, Zorth = Z; return; end

sz0 = size( Z );

sz = [ sz0( 1 ), sz0( 2 ) ];

if Nmodes > 1
   
    tilde_S = reshape( Z, [ sz( 1 ) * sz( 2 ), Nmodes ] );
    
    [ U, ~ ] = eig( tilde_S' * tilde_S );                       % U is unitary, i.e. U' = inv( U )
    
    S = tilde_S * U; 
%     S = tilde_S * ctranspose( U ); 
    
    Zorth = reshape( S, [ sz( 1 ), sz( 2 ), Nmodes ] );

end


% see "Breaking ambiguities in mixed state ptychography" paper

% see also "Reconstructing mode mixtures in the optical near-field" for fCDI ??


    
%{
    
[ u, s, v ] = svd( sDL.D );
% s = eye( sDL.Ni, sDL.Na );
% sDL.D = u * s * v';
sDL.D = u * v';
    
%}