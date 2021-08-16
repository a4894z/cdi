function [ TF ] = DMupdate_2Dsample_repmat( spos_conjP_exwv, ...
                                            spos_abs_P_abs2, ...
                                            phi,             ...
                                            TF,              ...
                                            probe,           ...
                                            ind,             ...
                                            rc,              ...
                                            Nspos )
             
% persistent spos_abs_P_abs2

%==============================
% vectorized ePIE sample update
%==============================
 
sum_conj_probe_phiM = reshape( sum( conj( probe ) .* phi, 3 ), [ rc, Nspos ] );     % sum over probe modes    

sum_abs2_probe = sum( abs( probe ) .^ 2, 3 );                                       % sum over probe modes    
sum_abs2_probe = repmat( sum_abs2_probe( : ), [ 1, Nspos ] );  
% sum_abs2_probe = sum_abs2_probe( : );  




% tic
% conj_probe = conj( probe );
% sum_conj_probe_phiM = reshape( sum( conj_probe .* phi, 3 ), [ rc, Nspos ] ); % sum over probe modes    
% 
% sum_abs2_probe = sum( probe .* conj_probe, 3 ); % sum over probe modes    
% sum_abs2_probe = repmat( sum_abs2_probe( : ), [ 1, Nspos ] );  
% % sum_abs2_probe = sum_abs2_probe( : );  
% toc

%========

%{

sz = length(TF);

spos_conjP_exwv = gpuArray.zeros( [ sz, Nspos ], 'single');
spos_abs_P_abs2 = gpuArray.zeros( [ sz, Nspos ], 'single');

%}


spos_conjP_exwv( ind ) = sum_conj_probe_phiM;     
spos_abs_P_abs2( ind ) = sum_abs2_probe;

spos_conjP_exwv = sum( spos_conjP_exwv - TF .* spos_abs_P_abs2, 2 );        % sum over scan positions
spos_abs_P_abs2 = sum( spos_abs_P_abs2, 2 );                                % sum over scan positions
TF = TF + spos_conjP_exwv / max( spos_abs_P_abs2 );  

% TF = TF + sum( spos_conjP_exwv - TF .* spos_abs_P_abs2, 2 ) / max( sum( spos_abs_P_abs2, 2 ) );   