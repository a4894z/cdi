function [ phi, va, lam, aDL ] = ADMMupdate_exwv_aDL_C_meas( Ca, lam, aDL, meas, sol )

%==================================================================================================
% ANALYSIS SPARSE COEFFICIENT UPDATE:

ua = Ca - lam;

% determine thresholding level for each column of X:
abs_ua = abs( ua );
temp1 = sort( abs_ua, 'descend' );
sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

% HARD SHRINKAGE
uaprime = ua .* ( abs_ua >= sel );  

%     % SOFT SHRINKAGE
%     abs_ua = abs( ua );
% %     ua = ( abs_ua - sel ) .* ( abs_ua >= sel ) .* ( ua ./ ( 1e-7 + abs_ua ));  
%     uaprime = ( abs_ua - sel ) .* ( abs_ua >= sel ) .* exp( 1i * angle( ua ));   

%==================================================================================================

va = uaprime + lam;

phi = trainingsetpatches2image( pinv( aDL.D ) * va, aDL, sol.sz );

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

[ Ip, aDL ] = image2trainingsetpatches( phi, aDL, sol.sz );

vaprime = aDL.D * Ip;

aDL.C = vaprime;

%==================================================================================================

lam = lam + 1.0 * ( uaprime - vaprime );
