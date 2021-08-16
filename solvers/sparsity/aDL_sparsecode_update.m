function [ Cs ] = aDL_sparsecode_update( R, aDL )

%==================================================================================================
% SPARSE ANALYSIS COEFFICIENT UPDATE

Cs = aDL.D * R;      % minimizer for min_Ca || Da * R - Ca ||^2_F

%==================================================================================================

% THRESHOLDING USING WHOLE SPARSE CODE MATRIX:

abs_C = abs( Cs );
thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * aDL.sparse_lvl );

Cs = Cs .* ( abs_C >= thresh );                                                             % hard thresh
% Cs = ( abs_C -  1.0 * thresh ) .* ( abs_C >= thresh ) .* exp( 1i * angle( Cs ));          % soft thresh

%==================================================================================================

% USING WHOLE SPARSE CODE MATRIX:


    
%     % determine thresholding level for each column of X:
%     abs_Ca = abs( Ca );
%     temp1 = sort( abs_Ca, 'descend' );
%     sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    
% 
%     %============================================
%     
% %     % hard thresholding
% %     Ca = Ca .* ( abs_Ca >= sel );  
%     
%     %============================================
%     
%     % soft thresholding
%     abs_Ca = abs( Ca );
%     Ca = ( abs_Ca -  1.0 * sel ) .* ( abs_Ca >= sel ) .* ( Ca ./ ( 1e-7 + abs_Ca ));  
%     
%     %============================================
% 
% %     % firm thresholding
% %     W2 = ( abs_Ca > sel );
% %     W1 = ( abs_Ca <= 0.5 * sel );
% %     W12 = not( W2 ) & not( W1 );
% % 
% %     Ca = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_Ca - 0.5 * sel ) ) .* ( Ca ./ ( 1e-7 + abs_Ca ) ) + Ca .* W2;
%   




%==================================================================================================