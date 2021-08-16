function [ phi, Cs, sDL ] = RAARupdate_exwvmeas_sDLdictionary( phi, Cs, sDL, meas, sol )


abs_C = abs( Cs );
thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * sDL.sparse_lvl );

Pi_s_Cs = Cs .* ( abs_C >= thresh );                                                             % hard thresh
% Pi_s_Cs = ( abs_C -  1.0 * thresh ) .* ( abs_C >= thresh ) .* exp( 1i * angle( Cs ));          % soft thresh







Ds = Ip * pinv( sDL.C );
% Ds = Ip * sDL.C';

% colsD_eq0 = ~any( Ds, 1 );
% Ds( :, colsD_eq0 ) = Ds0( :, colsD_eq0 );

if strcmp( sDL.is_orthog, 'orthog_yes' )
    
    [ u, ~, v ] = svd( Ds );
    Pi_s_Ds = u * v';

end


###














R_s_Cs = 2 * Pi_s_Cs - Cs;



 
 
 
 
%==================================================================================================

[ R_s_Cs_phi ] = trainingsetpatches2image( sDL.D * R_s_Cs, sDL, sol.sz.sz );

Pi_m_R_s_Cs_phi = enforce_2DTPAmeas( R_s_Cs_phi, meas, sol.measLPF, sol );

[ Pi_m_R_s_Cs_Ip, ~ ] = image2trainingsetpatches( Pi_m_R_s_Cs_phi, sDL, sol.sz.sz );

Pi_m_R_s_Cs = sDL.D' * Pi_m_R_s_Cs_Ip;                   % minimizer for min_Cs || Ds * Cs - Ip ||^2_F     !!!!!    s.t. Ds' Ds = I    !!!!!
% Pi_m_R_s_Cs = pinv( sDL.D ) * Pi_m_R_s_Cs_Ip;            % minimizer for min_Cs || Ds * Cs - Ip ||^2_F

R_m_R_s_Cs = 2 * Pi_m_R_s_Cs - R_s_Cs;

%==================================================================================================


[ phi ] = trainingsetpatches2image( sDL.D * Cs, sDL, sol.sz.sz );

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

[ Pi_m_Ip, ~ ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );

Pi_m_Cs = sDL.D' * Pi_m_Ip;                             % minimizer for min_Cs || Ds * Cs - Ip ||^2_F     !!!!!    s.t. Ds' Ds = I    !!!!!
% Pi_m_Cs = pinv( sDL.D ) * Pi_m_Ip;                      % minimizer for min_Cs || Ds * Cs - Ip ||^2_F

%==================================================================================================

bbeta = 0.5;

Cs =  bbeta * 0.5 * ( R_m_R_s_Cs + Cs ) + ( 1 - bbeta ) * Pi_m_Cs;



