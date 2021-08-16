function [ sample ] = DMupdate_2Dsample_loop_v2( sol )

tmp0 = zeros( sol.sample.sz.r, sol.sample.sz.c, 'single');
tmp1 = zeros( sol.sample.sz.r, sol.sample.sz.c, 'single');

sum_conj_probe_phiM = squeeze( sum( conj( sol.probe.P ) .* sol.phi, 3 )); 
sum_abs2_probe      = sum( abs( sol.probe.P ) .^ 2, 3 );

for ss = 1 : sol.spos.N % sol.spos.updateorder   % order ***DOESN'T*** matter here !!!
% for ss = sol.spos.updateorder   % order ***DOESN'T*** matter here !!!
    
    %==============
    
    rs = round( -1 * sol.spos.rs( ss, : ));

    %==============
    
    % TRY PAD ARRAY?
    shift_conjP_exwv = zeros( sol.sample.sz.r, sol.sample.sz.c, 'single');
    shift_abs_P_abs2 = zeros( sol.sample.sz.r, sol.sample.sz.c, 'single');
    
    shift_conjP_exwv( +1.0 * rs( 1 ) + sol.sample.vs.r, ...
                      +1.0 * rs( 2 ) + sol.sample.vs.c )      = sum_conj_probe_phiM( :, :, ss ); 

    shift_abs_P_abs2( +1.0 * rs( 1 ) + sol.sample.vs.r, ...
                      +1.0 * rs( 2 ) + sol.sample.vs.c )      = sum_abs2_probe;

    tmp0 = tmp0 + shift_conjP_exwv;
    tmp1 = tmp1 + shift_abs_P_abs2;
    
    %==============

end

% sample = sol.sample.TF + 10 * ( tmp0 - sol.sample.TF .* tmp1 ) / ( length(sol.spos.updateorder) * max( sum_abs2_probe( : )));       % 1e-7 vs eps wtf?

sample = sol.sample.TF + 1 * ( tmp0 - sol.sample.TF .* tmp1 ) / max( tmp1( : ));       % 1e-7 vs eps wtf?

% sample = ( tmp0 + sol.sample.TF .* tmp1 ) ./ ( 2 * ( 1e-7 + tmp1 ));       % 1e-7 vs eps wtf?














% 5;


% out.reverseStr = '';
% 
% for ss = sol.spos.updateorder
% 
%     out.msg = sprintf( 'DM sample update on scan position %d/%d', ss, sol.spos.N );
%     fprintf([out.reverseStr, out.msg]);
%     out.reverseStr = repmat(sprintf('\b'), 1, length(out.msg));
% 
% 
% end

