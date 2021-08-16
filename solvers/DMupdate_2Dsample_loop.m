function [ sample ] = DMupdate_2Dsample_loop( sol )

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

sample = tmp0 ./ ( tmp1 + 1e-7 );       % 1e-7 vs eps wtf?








% function [ sample ] = DMupdate_2Dsample_loop( phi, probe, updateorder, samsz, vs, spos_rs )
% 
% sum_conjP_exwv = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
% sum_abs_P_abs2 = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
% 
% shift_conjP_exwv = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
% shift_abs_P_abs2 = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
% 
% sum_conj_probe_phiM = squeeze( sum( conj( probe ) .* phi, 3 )); 
% sum_abs2_probe      = sum( abs( probe ) .^ 2, 3 );
% 
% 
% for ss = updateorder   % order ***DOESN'T*** matter here !!!
%     
%     %==============
%     
%     rs = round( -1 * spos_rs( ss, : ));
% 
%     %==============
% 
% %     shift_conjP_exwv = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
% %     shift_abs_P_abs2 = gpuArray( zeros( samsz( 1 ), samsz( 2 ), 'single'));
%     shift_conjP_exwv = 0 * shift_conjP_exwv;
%     shift_abs_P_abs2 = 0 * shift_abs_P_abs2;
%     
%     shift_conjP_exwv( +1.0 * rs( 1 ) + vs.r, ...
%                       +1.0 * rs( 2 ) + vs.c )      = sum_conj_probe_phiM( :, :, ss ); 
% 
%     shift_abs_P_abs2( +1.0 * rs( 1 ) + vs.r, ...
%                       +1.0 * rs( 2 ) + vs.c )      = sum_abs2_probe;
% 
%     sum_conjP_exwv = sum_conjP_exwv + shift_conjP_exwv;
%     sum_abs_P_abs2 = sum_abs_P_abs2 + shift_abs_P_abs2;
% 
%     %==============
% 
% %     sampleL = sum_conjP_exwv ./ ( 1e-8 + sum_abs_P_abs2 );
%     
% %     figure; imagesc( abs(sum_conjP_exwv))
% %     figure; imagesc( abs(sum_abs_P_abs2))
% %     figure; imagesc( abs(sampleL))
% 
% end
% 
% sample = sum_conjP_exwv ./ ( sum_abs_P_abs2 + eps );       % 1e-7 vs eps wtf?







































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

