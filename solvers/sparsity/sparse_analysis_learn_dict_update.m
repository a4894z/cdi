function [ Da, stepDa ] = sparse_analysis_learn_dict_update( Ca, Da, Ip, aDL )

stepDa = 0;

% R = zeros( size(Ca) );

tau = 5 / norm( Ca )^2;
%tau = 2 / norm( Ip )^2;
% tau = 1; %2 / norm( Ip )^2;


% for ii = 1 : 3
%     
%     %==========================================================================
%     %------------- minimizer ( wrt Da ) for  || Ca - Da * Ip ||_F^2 ------------
%     %==========================================================================
%     
%     
%     
% %     R = Da * Ip - Ca;
% %     Da = Da - tau * ( Da * Ip - Ca ) * Ip';
% 
% 
%     % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION
% 
%         
%         gradDa = ( Da * Ip - Ca ) * Ip';
% 
%         gamma_try = linspace( 0, tau, 11 );
% 
%         for ll = 1 : length( gamma_try )
% 
%             errDa( ll ) = norm( Ca - ( Da - gamma_try( ll ) * gradDa ) * Ip, 'fro' );
% 
%         end
% 
%         [ ~, I ] = min( errDa );
%         Da = Da - gamma_try( I ) * gradDa;
% 
%         figure; semilogy( gamma_try, errDa ); title( num2str( gamma_try( I ), 'step = %e' )); grid on
%         close all;
% 
%         stepDa = gamma_try( I );
% 
% 
% end
%     


    % for some reason, this doesn't work very well!?!?
    Da = Ca * pinv( Ip );
    
%     A*x = b
%     x1 = A\b 
%     x2 = pinv(A)*b
     
    %==============================================================================================
    
    
    
    


%     Da( isnan( Da )) = 0;
%     Da( isinf( Da )) = 0;
    
    
    
%     [U,S,V] = svd( Ip, 'econ' );
%     Sinv = 1 ./ S; 
%     Sinv( Sinv > 1e6 ) = 0;
%     Sinv2 = pinv(S);
%     D = Ca * ( V * pinv(S) * U' );






    %==============================================================================================
    
%     % this generally makes things worse
%     % make rows orthogonal:
%     T = Da.';
%     [ U, Do ] = eig( T' * T ); 
%     Da = ( T * U ).'; 

    % normalize dictionary atoms ( ROWS ) to be unit-norm. !!!!!!!!!!!!!!!!!!!!!!!!
    Da = Da ./ ( 1e-7 + repmat( sqrt( sum( abs( Da ) .^ 2, 2 )), [ 1, aDL.Ni ] ));  
    
%     % rescale so that each analysis dictionary row has max value of one
%     Da = Da ./ ( 1e-7 + repmat( max( abs( Da ) , [], 2 ), [ 1, aDL.Ni ] ));   

%     % rescale WHOLE MATRIX to be unit norm:
%     Da = Da / norm( Da, 'fro' );   

%     % rescale so that WHOLE MATRIX has max value of one
%     Da = Da / max( abs( Da(:) ) );   

    %==============================================================================================
    
    
    
    %{
    
    if any(isnan( Da(:) ))
        5; 
    end
    if any(isinf( Da(:) ))
        5; 
    end
    
    
    Da( isnan( Da )) = 0;
    Da( isinf( Da )) = 0;

    %}
    
    

    
