function [ phi, aDL ] = RAARupdate_exwv_aDL_C_meas( phi, aDL, meas, sol )
       
% 0.5 * bbeta * ( R_m[ R_s[ phi ]] + phi ) + ( 1 - bbeta ) * Pi_m[ phi ]

%==================================================================================================
% ORIGINAL IMAGE TO IMAGE PATCHES

[ Ip, aDL ] = image2trainingsetpatches( phi, aDL, sol.sz.sz );

%==================================================================================================
% SPARSE ANALYSIS COEFFICIENT UPDATE

aDL.C = aDL.D * Ip;      % minimizer for min_Ca || Da * Ip - Ca ||^2_F



% abs_C = abs( aDL.C );
% thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * aDL.sparse_lvl );
% 
% aDL.C = aDL.C .* ( abs_C >= thresh );                                                             % hard thresh
% % aDL.C = ( abs_C -  1.0 * thresh ) .* ( abs_C >= thresh ) .* exp( 1i * angle( aDL.C ));          % soft thresh
% 
% 
% 

%=======================

% determine thresholding level for each column of X:
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

%==================================================================================================

[ Pi_C_phi ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sol.sz.sz );

R_C_phi = 2 * Pi_C_phi - phi;

%==================================================================================================

Pi_m_R_C_phi = enforce_2DTPAmeas( R_C_phi, meas, sol.measLPF, sol );

R_m_R_C_phi = 2 * Pi_m_R_C_phi - R_C_phi;

%==================================================================================================

% Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
Pi_m_phi = enforce_2DTPAmeas( Pi_C_phi, meas, sol.measLPF, sol );

bbeta = 0.01;

phi = 0.5 * bbeta * ( R_m_R_C_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi;






