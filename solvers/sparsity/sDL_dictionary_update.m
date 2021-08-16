function [ Ds ] = sDL_dictionary_update( R, sDL )

%==================================================================================================
% minimizer for min_Ds || Ds * Cs - R ||^2_F           ( only applies if Ds is complete, not over/under )
% s.t. Ds' * Ds = identity matrix 

%%{

Ds0 = sDL.D;

Ds = R * pinv( sDL.C );
% Ds = R * sDL.C';

% colsD_eq0 = ~any( Ds, 1 );
% Ds( :, colsD_eq0 ) = Ds0( :, colsD_eq0 );

% Ds = real( Ds );

if strcmp( sDL.is_orthog, 'orthog_yes' )
    
    [ u, ~, v ] = svd( Ds );
    Ds = u * v';

else

%     Ds = Ds ./ ( 1e-7 + repmat( max( abs( Ds )), [ size( Ds, 1 ), 1 ] ));                 % rescale so that each synthesis dictionary component has max value of one

%     tmp0 = abs( Ds ) > 1; Ds( tmp0 ) = 1 * exp( 1i * angle( Ds( tmp0 )));                 % project so that each abs() synthesis dictionary component has max value of one

    Ds = Ds ./ ( 1e-7 + repmat( sqrt( sum( abs( Ds ) .^ 2, 1 )), [ size( Ds, 1 ), 1 ] ));   % normalize dictionary atoms ( COLUMNS ) to be unit-norm. 

%     [ Ds ] = scale_matrixcolumns_phase( Ds, [ 0, 1 ] );
%     [ Ds ] = scale_matrixcolumns_magnitude( Ds, 1.0 * [ -pi, pi ] );
% 
end

% [ u, s, v ] = svd( Ds );
% s = eye( size( Ds, 1 ), size( Ds, 2 ) );
% Ds = u * s * v';

%}

%==================================================================================================
% minimizer for min_Ds || Ds * Cs - R ||^2_F,      s.t. cols( abs( Ds )) <= 1, or ???
% 
% Ds = R * pinv( Cs );    
% % Ds = R * Cs';        % minimizer for min_Ds || Ds * Cs - R ||^2_F           !!!! ASSUMING !!!! Cs' = pinv( Cs )
% 
% % colsD_eq0 = ~any( Ds, 1 );
% % Ds( :, colsD_eq0 ) = Ds0( :, colsD_eq0 );

%=======================

%==================================================================================================

% Ds( isnan( Ds )) = 0;
% Ds( isinf( Ds )) = 0;



%{

% WHAT TO DO ABOUT ZERO COLUMNS IN DICTIONARY ???
% EFFECTIVELY MAKES IT UNDERCOMPLETE ??

Ds0 = Ds;
colsD_eq0 = ~any( Ds0, 1 );
Ds0( :, colsD_eq0 ) = []; 

Cs0 = Cs;
Cs0( colsD_eq0, :  ) = []; 

%}


%==================================================================================================

% return
% 
% 
% 
% stepDs = 1;
% 
% % R = zeros( size( Cs ) );
% 
% tau = 3 / norm( Cs )^2;
% %tau = 1; 
% 
% CsH = Cs';
% 
% z = 0;
% 
% for ii = 1 : 3 %sDL.ndict
%     
%     %==========================================================================
%     %------------- minimizer ( wrt Ds ) for  || Ds * Cs - R ||_F^2 ------------
%     %==========================================================================
%     
%     
% %     % projected steepest descent:
% %     Ds = Ds - tau * ( Ds * Cs - R ) * CsH;
% 
% 
% 
% 
% 
% 
% 
%     % TRY EXACT LINE SEARCH USING THIS GRADIENT AND ALSO WITH L0 OR L1 REGULARIZATION
%     
%     gradDs = ( Ds * Cs - R ) * CsH;
%     
%     gamma_try = linspace( 0, tau, 11 );
%     
%     for ll = 1 : length( gamma_try )
%         
%         errDs( ll ) = norm( ( Ds - gamma_try( ll ) * gradDs ) * Cs - R, 'fro' );
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

