function [ TF ] = mrPIEupdate_sample( probe, TF, phiM, sol, expt )
  
%==================================================================================================

% sum_abs_probe_abs2 = zeros( sol.sz.r, sol.sz.c, 'single');
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     sum_abs_probe_abs2 = sum_abs_probe_abs2 + abs( probe( :, :, pp )) .^ 2; 
% 
% end

sum_abs_probe_abs2 = sum( abs( probe ) .^ 2, 3 );

max_sum_abs_probe_abs2 = max( sum_abs_probe_abs2( : ));

%==================================================================================================

if ~isfield( sol.spos, 'updateorder' )
    
    updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));
    
else
    
    updateorder = sol.spos.updateorder;
    
end

%==================================================================================================

conj_probe = conj( probe );

rPIE_alpha = 0.90;
update_term_T_denom = ( rPIE_alpha * max_sum_abs_probe_abs2 + ( 1 - rPIE_alpha ) * sum_abs_probe_abs2 );

tt = 0;
Tmom = 10;
eta_a = 0.7;
eta_b = 0.7;
gamma_m = 0.3;
TF_j_plus_1_minusT = TF;
nu_j_minus_T = zeros( sol.sample.sz.sz, 'single');

for ss = updateorder       % order * DOES * matter here !!!
  
    %========================================
        
% PRINT OUT % DONE, EVERY 10%

    %========================================

    rs = +sol.spos.rs( ss, : );                 % !!!!!!!!!!!
    
    [ TFview ] = getview_2DsampleTF( TF, sol.sample.vs, rs, sol.spos.shifttype );

    update_term_T = sum( conj_probe .* ( phiM( :, :, :, ss ) - probe .* TFview ), 3 );

%     TF( round( -1.0 * rs( 1 ) + sol.sample.vs( 1, : ) ), ...
%         round( -1.0 * rs( 2 ) + sol.sample.vs( 2, : ) )) = TFview + gamma_m * update_term_T / max_abs_probe_abs2;  % !!!!!!
    
    update_term_T = update_term_T ./ update_term_T_denom;
    
    TF( round( -1.0 * rs( 1 ) + sol.sample.vs.r ), ...
        round( -1.0 * rs( 2 ) + sol.sample.vs.c )) = TFview + gamma_m * update_term_T;
    
    %========================================
    
    TF_prime = TF;
    tt = tt + 1;

    if ( mod( tt, Tmom ) == 0 )
        
        nu_j = eta_a * nu_j_minus_T + TF_prime - TF_j_plus_1_minusT;

%         TF = TF + nu;
        TF = TF_prime + eta_b * nu_j;
         
        tt = 0;
        nu_j_minus_T = nu_j;
        TF_j_plus_1_minusT = TF;
        
    end
   
end
