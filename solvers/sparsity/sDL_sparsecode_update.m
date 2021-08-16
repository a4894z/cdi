function [ Cs ] = sDL_sparsecode_update( R, sDL )

%==================================================================================================

if strcmp( sDL.is_orthog, 'orthog_yes' )
    
    Cs = sDL.D' * R;                   % minimizer for min_Cs || Ds * Cs - R ||^2_F     !!!!!    s.t. Ds' Ds = I    !!!!!
       
else
    
    Cs = pinv( sDL.D ) * R;            % minimizer for min_Cs || Ds * Cs - R ||^2_F
    
end

% norm( sDL.D * pinv( sDL.D ) * R - R, 'fro' )

%==================================================================================================
%---------------------- ENFORCE SPARSITY LEVEL ON IMAGE PATCH SPARSE CODES ------------------------
%==================================================================================================

abs_C = abs( Cs );
thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * sDL.sparse_lvl );

Cs = Cs .* ( abs_C >= thresh );                                                             % hard thresh
% Cs = ( abs_C -  1.0 * thresh ) .* ( abs_C >= thresh ) .* exp( 1i * angle( Cs ));          % soft thresh

%=======================

% re_Cs = real( Cs );
% im_Cs = imag( Cs );
% thresh_re = find_thresh_from_sparsitylevel( re_Cs, numel( re_Cs ) * sDL.sparse_lvl );
% thresh_im = find_thresh_from_sparsitylevel( im_Cs, numel( im_Cs ) * sDL.sparse_lvl );
% 
% re_Cs = re_Cs .* ( re_Cs >= thresh_re );  
% im_Cs = im_Cs .* ( im_Cs >= thresh_im );  
% Cs = re_Cs + 1i * im_Cs;
% 
% % re_Cs = ( re_Cs -  1.0 * thresh_re ) .* ( re_Cs >= thresh_re ) .* sign( re_Cs );  
% % im_Cs = ( im_Cs -  1.0 * thresh_im ) .* ( im_Cs >= thresh_im ) .* sign( im_Cs );  
% % Cs = re_Cs + 1i * im_Cs;

%==================================================================================================
% determine thresholding level for each column of sparse code matrix:


% abs_C = abs( Cs );
% temp1 = sort( abs_C, 'descend' );
% sel = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );  
% 
% % hard thresholding
% Cs = Cs .* ( abs_C >= sel );  

% % soft thresholding
% % Cs = ( abs_C -  1.0 * sel ) .* ( abs_C >= sel ) .* ( Cs ./ ( 1e-7 + abs_C ));  
% Cs = ( abs_C -  1.0 * sel ) .* ( abs_C >= sel ) .* exp( 1i * angle( Cs ));  
% 
% 
% % firm thresholding
% W2 = ( abs_C > sel );
% W1 = ( abs_C <= 0.5 * sel );
% W12 = not( W2 ) & not( W1 );
% Cs = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_C - 0.5 * sel ) ) .* ( Cs ./ ( 1e-7 + abs_C )) + Cs .* W2;




%{







% % sel = 0.005 + 0 * sel;
% 
% % perform thresholding:
% if strcmp( sDL.thresh, 'hard' )
%     
%     % hard thresholding
%     Cs = Cs .* ( abs_C >= sel );  
% 
% elseif strcmp( sDL.thresh, 'soft' )
% 
%     % soft thresholding
% %     Cs = ( abs_C -  1.0 * sel ) .* ( abs_C >= sel ) .* ( Cs ./ ( 1e-7 + abs_C ));  
%     Cs = ( abs_C -  1.0 * sel ) .* ( abs_C >= sel ) .* exp( 1i * angle( Cs ));  
%     
% elseif strcmp( sDL.thresh, 'firm' )
% 
%     % firm thresholding
%     W2 = ( abs_C > sel );
%     W1 = ( abs_C <= 0.5 * sel );
%     W12 = not( W2 ) & not( W1 );
% 
%     Cs = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_C - 0.5 * sel ) ) .* ( Cs ./ ( 1e-7 + abs_C )) + Cs .* W2;
% 
% else
%     
%     error('SPECIFY PROPER THRESHOLDING METHOD, EXITING.')
%     
% end

%}

%==================================================================================================

%{

% determine thresholding level for each column of sparse code matrix:
re_C = real( Cs );
temp1 = sort( re_C, 'descend' );
sel_re = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );  

im_C = imag( Cs );
temp1 = sort( im_C, 'descend' );
sel_im = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

% perform thresholding:
if strcmp( sDL.thresh, 'hard' )
    
    % hard thresholding
    re_C = re_C .* ( re_C >= sel_re );  
    im_C = im_C .* ( im_C >= sel_im );  
     
elseif strcmp( sDL.thresh, 'soft' )

    % soft thresholding
    re_C = ( re_C -  1.0 * sel_re ) .* ( re_C >= sel_re ) .* sign( re_C );  
    im_C = ( im_C -  1.0 * sel_im ) .* ( im_C >= sel_im ) .* sign( im_C );  
    
elseif strcmp( sDL.thresh, 'firm' )

    % firm thresholding
    W2 = ( re_C > sel_re );
    W1 = ( re_C <= 0.5 * sel_re );
    W12 = not( W2 ) & not( W1 );
    re_C = W12 .* ( ( sel_re ./ ( sel_re - 0.5 * sel_re ) ) .* ( re_C - 0.5 * sel_re ) ) .* sign( re_C ) + re_C .* W2;

    W2 = ( im_C > sel_im );
    W1 = ( im_C <= 0.5 * sel_im );
    W12 = not( W2 ) & not( W1 );
    im_C = W12 .* ( ( sel_im ./ ( sel_im - 0.5 * sel_im ) ) .* ( im_C - 0.5 * sel_im ) ) .* sign( im_C ) + im_C .* W2;
    
else
    
    error('SPECIFY PROPER THRESHOLDING METHOD, EXITING.')
    
end

Cs = re_C + 1i * im_C;

%}

%==================================================================================================


% Cs = Cs ./ ( 1e-7 + repmat( max( abs( Cs ) , [], 1 ), [ size( Cs, 1 ), 1 ] ));            % rescale so that each analysis sparse code COLUMN has max abs value of one

% max_abs = 1; tmp0 = ( abs( Cs ) > max_abs );
% Cs( tmp0 ) = max_abs * exp( 1i * angle( Cs( tmp0 )));                                     % enforce max sparse code value

%==================================================================================================
% steepest descent with line search:

% DsH = sDL.D';
% gradCs = DsH * ( sDL.D * sDL.Cs - R );
% 
% ggamma = norm( Cs, 'fro' ) / norm( gradCs, 'fro' );
% Cs = sDL.Cs - 0.1 * ggamma * gradCs;

%{

stepCs = 1;

gamma = 3 / norm( sDL.D )^2;  % norm( sDL.D ) ~ max(svd( sDL.D ))
%gamma = 1; 

DsH = sDL.D';

z = 0;



    
    for ii = 1 : 1




    %     % projected steepest descent:
    %     R = sDL.D * Cs - R;
    %     Cs = Cs - gamma * DsH * R;
    % 
    % %     stepCs = ( norm( Cs, 'fro' ) / norm( DsH * R, 'fro' ));
    % %     Cs = Cs - stepCs * DsH * R;









        % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION

        gradCs = DsH * ( sDL.D * sDL.Cs - R );
        %gamma = 3 / norm( gradCs, 'fro' )^1;

        gamma_try = linspace( 0, gamma, 21 );

        for ll = 1 : length( gamma_try )

            errCs( ll ) = norm( sDL.D * ( sDL.Cs - gamma_try( ll ) * gradCs ) - R, 'fro' );

        end

        [ ~, I ] = min( errCs );

        Cs = sDL.Cs - 1 * gamma_try( I ) * gradCs;

        stepCs = gamma_try( I );  
        5;


        figure; 
        semilogy( gamma_try, errCs ); title( num2str( gamma_try( I ), 'step = %e' )); grid on 
        close all;

    end

%}

%==================================================================================================



