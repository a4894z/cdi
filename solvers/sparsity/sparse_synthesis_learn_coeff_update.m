function [ Cs, stepCs ] = sparse_synthesis_learn_coeff_update( Cs, Ds, Ip, sDL )

% R = zeros( size( Cs ));

stepCs = 1;

gamma = 3 / norm( Ds )^2;  % norm( Ds ) ~ max(svd( Ds ))
%gamma = 1; 

DsH = Ds';

z = 0;


    %==========================================================================
    %------------- minimizer ( wrt Cs ) for  || Ds * Cs - Ip ||_F^2 ------------
    %==========================================================================

    
    for ii = 1 : 3 %sDL.ncoef




    %     % projected steepest descent:
    %     R = Ds * Cs - Ip;
    %     Cs = Cs - gamma * DsH * R;
    % 
    % %     stepCs = ( norm( Cs, 'fro' ) / norm( DsH * R, 'fro' ));
    % %     Cs = Cs - stepCs * DsH * R;









        % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION

        gradCs = DsH * ( Ds * Cs - Ip );
        %gamma = 3 / norm( gradCs, 'fro' )^1;

        gamma_try = linspace( 0, gamma, 11 );

        for ll = 1 : length( gamma_try )

            errCs( ll ) = norm( Ds * ( Cs - gamma_try( ll ) * gradCs ) - Ip, 'fro' );

        end

        [ ~, I ] = min( errCs );

        Cs = Cs - 1 * gamma_try( I ) * gradCs;

        stepCs = gamma_try( I );  
        5;


    %     figure; 
    %     semilogy( gamma_try, errCs ); title( num2str( gamma_try( I ), 'step = %e' )); grid on 
    %     close all;

    end


    %===========================

    % for some reason, this doesn't work very well!?!?

%     tol = max( size( Ds )) * eps( norm( Ds ));
%     tol = 1e-2;
%     Cs = pinv( Ds, tol ) * Ip;
    Cs = pinv( Ds ) * Ip;

%     A*x = b
%     x1 = A\b 
%     x2 = pinv(A)*b

    %==============================================================================================
    
    % determine thresholding level for each column of Cs:
    abs_Cs = abs( Cs );
    temp1 = sort( abs_Cs, 'descend' );
    sel = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

    %============================================
    
%     % hard thresholding
%     Cs = Cs .* ( abs_Cs >= sel );  
    
    %============================================
    
    % soft thresholding
    abs_Cs = abs( Cs );
    %Cs = ( abs_Cs - sel ) .* ( abs( Cs ) >= sel ) .* exp( 1i * angle( Cs ));  
    Cs = ( abs_Cs - 1 * sel ) .* ( abs( Cs ) >= sel ) .* ( Cs ./ ( 1e-7 + abs_Cs ) );  

    %============================================

%     % firm thresholding
%     W2 = ( abs_Cs > sel );
%     W1 = ( abs_Cs <= 0.5 * sel );
%     W12 = not( W2 ) & not( W1 );
% 
%     Cs = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_Cs - 0.5 * sel ) ) .* ( Cs ./ ( 1e-7 + abs_Cs ) ) + Cs .* W2;

    %==============================================================================================
    
%     % rescale so that each synthesis sparse code COLUMN has max abs value of one
%     Cs = Cs ./ ( 1e-7 + repmat( max( abs( Cs )), [ sDL.Na, 1 ] ));
    
    % ENFORCE MAX SPARSE CODE VALUE, IF BIGGER THAN 1, PROJECT ONTO 1
    tmp0 = abs( Cs ) > 1;
    Cs( tmp0 ) = 1 * exp( 1i * angle( Cs( tmp0 )));
    
    %==============================================================================================
    
    %{
    
    if any(isnan( Cs(:) ))
        5; 
    end
    if any(isinf( Cs(:) ))
        5; 
    end


    Cs( isnan( Cs )) = 0;
    Cs( isinf( Cs )) = 0;
    
    %}






