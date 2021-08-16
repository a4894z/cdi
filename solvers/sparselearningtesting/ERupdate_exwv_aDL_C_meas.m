function [ phi, aDL ] = ERupdate_exwv_aDL_C_meas( phi, aDL, meas, sol )

%==================================================================================================
% EXIT WAVE MEASUREMENT CONSTRAINT UPDATE:

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

% tmp0 = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% phi = 2 * tmp0 - phi;

% phi = phi .* sol.probe.S;

%==================================================================================================
% SPARSE ANALYSIS COEFFICIENT UPDATE

[ Ip, aDL ] = image2trainingsetpatches( phi, aDL, sol.sz );

%=======================

aDL.C = aDL.D * Ip;      % minimizer for min_Ca || Da * Ip - Ca ||^2_F

%=======================

% determine thresholding level for each column of sparse code matrix:
abs_C = abs( aDL.C );
temp1 = sort( abs_C, 'descend' );
sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

%=======================

% hard thresholding
aDL.C = aDL.C .* ( abs_C >= sel );  

%=======================

% % soft thresholding
% aDL.C = ( abs_C -  1.0 * sel ) .* ( abs_C >= sel ) .* ( aDL.C ./ ( 1e-7 + abs_C ));  

%=======================

% % firm thresholding
% W2 = ( abs_C > sel );
% W1 = ( abs_C <= 0.5 * sel );
% W12 = not( W2 ) & not( W1 );
% 
% aDL.C = W12 .* ( ( sel ./ ( sel - 0.5 * sel ) ) .* ( abs_C - 0.5 * sel ) ) .* ( aDL.C ./ ( 1e-7 + abs_C )) + aDL.C .* W2;

%=======================

% % rescale so that each analysis sparse code COLUMN has max abs value of one
% aDL.C = aDL.C ./ ( 1e-7 + repmat( max( abs( aDL.C ) , [], 1 ), [ aDL.Na, 1 ] ));   

%     % enforce max sparse code value
%     tmp0 = abs( Ca ) > 1;
%     aDL.C( tmp0 ) = 1 * exp( 1i * angle( aDL.C( tmp0 )));

%=======================

[ phi ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sol.sz );

%==================================================================================================


