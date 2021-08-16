function [ TF ] = ePIEupdate_sampleOLD( probe, TF, phiM, sol, expt )

%==================================================================================================
    
sum_abs_probe_abs2 = zeros( sol.sz.r, sol.sz.c, 'single');

for pp = 1 : sol.probe.scpm.N
    
    %pad_probe( :, :, pp ) = padarray( probe( :, :, pp ), pad_amount ); 
    sum_abs_probe_abs2 = sum_abs_probe_abs2 + abs( probe( :, :, pp )) .^ 2; 

end

max_abs_probe_abs2 = max( sum_abs_probe_abs2( : ));

%==================================================================================================

pad_amount = double( round( [ 0.5 * ( sol.sample.sz.r - sol.sz.r ), 0.5 * ( sol.sample.sz.c - sol.sz.c ) ] ));

% pad_sum_abs_probe_abs2 = padarray( sum_abs_probe_abs2, pad_amount );
% pad_sum_abs_probe_abs2 = zeropadarray( sum_abs_probe_abs2, pad_amount );

% out.reverseStr = '';

% spos_updateorder = randperm( length( sol.spos.rs ));

for ss = sol.spos.updateorder
    
    rs = -1 * sol.spos.rs( ss, : );                 % !!!!!!!!!!!

    %========================================

%     out.msg = sprintf( 'ePIE sample update on scan position %d/%d', ss, sol.spos.N );
%     fprintf( [ out.reverseStr, out.msg ] );
%     out.reverseStr = repmat( sprintf( '\b' ), 1, length( out.msg ));

    %========================================
        
%     spos = -1 * sol.spos.rs( ss , : );
    update_term_T = zeros( sol.sample.sz.r, sol.sample.sz.c, 'single');
    
    for pp = 1 : sol.probe.scpm.N
        
        pad_probe = padarray( probe( :, :, pp ), pad_amount );          % can store this outside of spos loop, don't need to shift and zero pad everytime
%         pad_probe = zeropadarray( probe( :, :, pp ), pad_amount );          
        pad_shft_probe = circshift( pad_probe, round( rs ));
%         pad_shft_probe = circshift( pad_probe, rs );
        
        pad_phiM = padarray( phiM( :, :, pp, ss ), pad_amount ); 
%         pad_phiM = zeropadarray( phiM( :, :, pp, ss ), pad_amount );
        pad_shft_phiM  = circshift( pad_phiM, round( rs ));

        update_term_T = update_term_T + conj( pad_shft_probe ) .* ( pad_shft_phiM - pad_shft_probe .* TF );

    end

%     TF = TF + 1.0 * update_term_T ./ ( 0.5 * max_abs_probe_abs2 + 0.5 * pad_sum_abs_probe_abs2 );
    TF = TF + 1.0 * update_term_T / max_abs_probe_abs2;
    
end
