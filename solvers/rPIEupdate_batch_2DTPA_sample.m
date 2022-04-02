function [ T ] = rPIEupdate_batch_2DTPA_sample( psi,   ...
                                                T,     ...
                                                probe, ...
                                                ind,   ...
                                                rc,    ...
                                                Nspos, ...
                                                rPIE_alpha )
                                     
                                     
                                       
%==============================
% vectorized rPIE sample update
%==============================
 
sum_conj_probe_psiM = reshape( sum( conj( probe ) .* psi, 3 ), [ rc, Nspos ] ); % sum over probe modes    

sum_abs2_probe = sum( abs( probe ) .^ 2, 3 );                                   % sum over probe modes    
sum_abs2_probe = repmat( sum_abs2_probe( : ), [ 1, Nspos ] );  

%========

sz = length( T );

spos_conjP_exwv = gpuArray.zeros( [ sz, Nspos ], 'single');
spos_abs_P_abs2 = gpuArray.zeros( [ sz, Nspos ], 'single');

spos_conjP_exwv( ind ) = sum_conj_probe_psiM;     
spos_abs_P_abs2( ind ) = sum_abs2_probe;

spos_conjP_exwv = sum( spos_conjP_exwv - T .* spos_abs_P_abs2, 2 );        % sum over scan positions
spos_abs_P_abs2 = sum( spos_abs_P_abs2, 2 );                                % sum over scan positions

T = T + spos_conjP_exwv ./ ( rPIE_alpha * max( spos_abs_P_abs2 ) + ( 1 - rPIE_alpha ) * spos_abs_P_abs2 ); 







% figure; 
% imagesc( (max( max( sum_abs2_probe )) - sum_abs2_probe ))
% daspect([1 1 1])
% colormap bone
% 
% 

% figure; 
% imagesc( (max( max( reshape( spos_abs_P_abs2, [1280, 1280] ))) - reshape( spos_abs_P_abs2, [1280, 1280] )))
% daspect([1 1 1])
% colormap bone