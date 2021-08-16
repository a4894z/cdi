function [ TF ] = rPIEupdate_sample( probe, TF, phiM, sol, expt )
  
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

rPIE_alpha = 0.1;

update_term_T_denom = ( rPIE_alpha * max_sum_abs_probe_abs2 + ( 1 - rPIE_alpha ) * sum_abs_probe_abs2 );

for ss = updateorder       % order * DOES * matter here !!!
  
    %========================================
        
% PRINT OUT % DONE, EVERY 10%

    %========================================

    rs = +sol.spos.rs( ss, : );                 % !!!!!!!!!!!
    
    [ TFview ] = getview_2DsampleTF( TF, sol.sample.vs, rs, sol.spos.shifttype );

    update_term_T = sum( conj_probe .* ( phiM( :, :, :, ss ) - probe .* TFview ), 3 );
    
%     update_term_T = update_term_T ./ ( rPIE_alpha * max_abs_probe_abs2 + one_minus_rPIE_alpha * sum_abs_probe_abs2 );
    update_term_T = update_term_T ./ update_term_T_denom;
    
    TF( round( -1.0 * rs( 1 ) + sol.sample.vs.r ), ...
        round( -1.0 * rs( 2 ) + sol.sample.vs.c )) = TFview + update_term_T;
    
    %========================================
    
end
