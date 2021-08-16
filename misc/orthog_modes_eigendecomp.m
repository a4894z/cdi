% function [ Zorth ] = orthog_modes_eigendecomp( Z, sz, Nmodes )
function [ Zorth ] = orthog_modes_eigendecomp( Z )

Nmodes = size( Z, 3 );
if Nmodes == 1, Zorth = Z; return; end

sz0 = size( Z );
sz = [ sz0( 1 ), sz0( 2 ) ];

if Nmodes > 1
   
    tilde_S = reshape( Z, [ sz( 1 ) * sz( 2 ), Nmodes ] );
    
    %[ U, D, V ] = eig( tilde_S' * tilde_S ); 
    [ U, ~ ] = eig( tilde_S' * tilde_S ); 
%     [ U, ~ ] = eig( tilde_S' * tilde_S, 'nobalance' );  % 'chol', 'qz'
    
    S = tilde_S * U; 

    Zorth = reshape( S, [ sz( 1 ), sz( 2 ), Nmodes ] );

    
%{
    
[ u, s, v ] = svd( sDL.D );
% s = eye( sDL.Ni, sDL.Na );
% sDL.D = u * s * v';
sDL.D = u * v';
    
%}
    
end


% see "Breaking ambiguities in mixed state ptychography" paper

% see also "Reconstructing mode mixtures in the optical near-field" for fCDI ??

