function [ Ds, stepDs ] = sparse_synthesis_learn_dict_update( Cs, Ds, Ip, sDL )

stepDs = 1;

% R = zeros( size( Cs ) );

% tau = 3 / norm( Cs )^2;
% %tau = 1; 
% 
% CsH = Cs';
% 
% z = 0;


% for ii = 1 : 3 %sDL.ndict
%     
%     %==========================================================================
%     %------------- minimizer ( wrt Ds ) for  || Ds * Cs - Ip ||_F^2 ------------
%     %==========================================================================
%     
%     
% %     % projected steepest descent:
% %     Ds = Ds - tau * ( Ds * Cs - Ip ) * CsH;
% 
% 
% 
% 
% 
% 
% 
%     % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION
%     
%     gradDs = ( Ds * Cs - Ip ) * CsH;
%     
%     gamma_try = linspace( 0, tau, 11 );
%     
%     for ll = 1 : length( gamma_try )
%         
%         errDs( ll ) = norm( ( Ds - gamma_try( ll ) * gradDs ) * Cs - Ip, 'fro' );
%         
%     end
%     
%     [ ~, I ] = min( errDs );
%     
% %     z = 0.999 * z + gradDs;
% %     Ds = Ds - gamma_try( I ) * z;
%     
%     Ds = Ds - 1 * gamma_try( I ) * gradDs;
%     
%     
%     stepDs = gamma_try( I );
%     
%     
%     
% %     figure; 
% %     semilogy( gamma_try, errDs ); title( num2str( gamma_try( I ), 'step = %e' )); grid on
% %     close all;
% 
%     
%     % TRY MOMENTUM METHODS
% 
% 
% 
% 
%     
% 
% end

    
    %===========================


    % for some reason, this doesn't work very well!?!?
    Ds = Ip * pinv( Cs );

%     A*x = b
%     x1 = A\b 
%     x2 = pinv(A)*b




    
%     [ U, S, V ] = svd( Cs );
%     Ds = Ip * ( V * pinv( S ) * U' );



    %==============================================================================================
    
%     % normalize dictionary atoms ( COLUMNS ) to be unit-norm. 
%     Ds = Ds ./ ( 1e-7 + repmat( sqrt( sum( abs( Ds ) .^ 2, 1 )), [ sDL.Ni, 1 ] ));  
%     Ds = Ds ./ ( 1e-7 + repmat( norm( Ds, 'fro' ), [ sDL.Ni, 1 ] ));  
    
%     % rescale so that each synthesis dictionary component has max value of one
%     Ds = Ds ./ ( 1e-7 + repmat( max( abs( Ds )), [ sDL.Ni, 1 ] ));

    %==============================================================================================
    
    %{
    
    if any(isnan( Ds(:) ))
        5; 
    end
    if any(isinf( Ds(:) ))
        5; 
    end

    Ds( isnan( Ds )) = 0;
    Ds( isinf( Ds )) = 0;
    
    %}
    

%     % orthogonal/orthonormal dictionary columns:
%     [ U, E ] = eig( Ds' * Ds ); 
%     %Ds = Ds * U ./ ( 1e-6 + repmat( sqrt( diag( E )).' , [ sDL.Ni, 1 ] ) ); 
%     Ds = Ds * U; 
%     

 
%{

% gram-schmidt orthogonalization

[~, Do1] = mgsog( Ds );
[~, Do2] = gsog( Ds );
[~, Do3] = gson( Ds );
[~, Do4] = mgsog( Ds );
[~, Do4] = mgsog( Ds );
[ Q, R ] = qr( Ds );

%     Ds = orth( Ds );
    

%}

    
    



