function [ Ca, stepCa ] = sparse_analysis_learn_coeff_update( Ca, Da, Ip, aDL )

% gamma = 2 / norm( Da )^2;  % norm( D ) ~ max(svd( D )) if D is a matrix (max singular value)

% R = zeros( size( Ca ));
tau = 3; %2 / norm( Ip )^2;
stepCa = 0;

for ii = 1 : aDL.ncoef

    %==========================================================================
    %------------- minimizer ( wrt Ca ) for  || Ca - Da * Ip ||_F^2 ------------
    %==========================================================================

%     % for some reason, this doesn't work very well!?!?
%     R = Ca - Da * Ip;
%     temp0 = Ca - gamma * R;





%     % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION
%     
%     gradCa = Ca - Da * Ip;
%     gradCa = Ca - Da * Ip + ttau * Ca ./ abs( Ca ) ;
%     
%     gamma_try = linspace( -tau, tau, 101 );
%     
%     for ll = 1 : length( gamma_try )
%         
%         errCa( ll ) = norm( ( Ca - gamma_try( ll ) * gradCa ) - Da * Ip, 'fro' );
% %         errCa( ll ) = norm( ( Ca - gamma_try( ll ) * gradCa ) - Da * Ip, 'fro' ) + sum( sum( abs( ( Ca - gamma_try( ll ) * gradCa ) )));
%     end
%     
%     [ ~, I ] = min( errCa );
%     Ca = Ca - gamma_try( I ) * gradCa;
%     
%     %{
%     
%     figure; 
%     semilogy( gamma_try, errCa ); 
%     title( num2str( gamma_try( I ), 'step = %e' )); 
%     grid on
%     
%     close all;
% 
%     %}
%     
%     stepCa = gamma_try( I );







   
    Ca = Da * Ip;      % minimizer for min_Ca || Da * Ip - Ca ||^2_F
   
    %==============================================================================================
    
    % determine thresholding level for each column of X:
    abs_Ca = abs( Ca );
    temp1 = sort( abs_Ca, 'descend' );
    sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

    %============================================
    
%     % hard thresholding
%     Ca = Ca .* ( abs_Ca >= sel );  
    
    %============================================
    
    % soft thresholding
    abs_Ca = abs( Ca );
    Ca = ( abs_Ca -  1.0 * sel ) .* ( abs_Ca >= sel ) .* ( Ca ./ ( 1e-7 + abs_Ca ));  
    
    %============================================

%     % firm thresholding
%     W2 = ( abs_Ca > sel );
%     W1 = ( abs_Ca <= 0.5 * sel );
%     W12 = not( W2 ) & not( W1 );
% 
%     Ca = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_Ca - 0.5 * sel ) ) .* ( Ca ./ ( 1e-7 + abs_Ca ) ) + Ca .* W2;
  
    %==============================================================================================
   
%     % rescale so that each analysis sparse code COLUMN has max abs value of one
%     Ca = Ca ./ ( 1e-7 + repmat( max( abs( Ca ) , [], 1 ), [ aDL.p, 1 ] ));   

%     % enforce max sparse code value
%     tmp0 = abs( Ca ) > 1;
%     Ca( tmp0 ) = 1 * exp( 1i * angle( Ca( tmp0 )));
    
    %==============================================================================================
    
    
end

