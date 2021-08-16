function [ P ] = DMupdate_probemodes( T, phi, sol )

% tmp99          = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );
sum_conjT_exwv = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );
sum_abs_T_abs2 = zeros( sol.sz.r, sol.sz.c, 'single' );

% out.reverseStr = '';

%==================================================================================================

if ~isfield( sol.spos, 'updateorder' )
    
    updateorder = randperm( length( sol.spos.indxsubset ), 1.0 * length( sol.spos.indxsubset ));
else
    
    updateorder = sol.spos.updateorder;
    
end

%==================================================================================================
% tic
for ss = updateorder           % order * DOESN'T * matter here !!!

%     rs = +1 * sol.spos.rs( ss, : );    % !!!!!!!!!!
%     rs = -1 * sol.spos.rs( ss, : );      
    
    %============================================

%     out.msg = sprintf( 'DM probe update on scan position %d/%d', ss, sol.spos.N );
%     fprintf([out.reverseStr, out.msg]);
%     out.reverseStr = repmat(sprintf('\b'), 1, length(out.msg));

    %============================================
  
    TFview = getview_2DsampleTF( T, sol.sample.vs, sol.spos.rs( ss, : ), sol.spos.shifttype );
    
    sum_abs_T_abs2 = sum_abs_T_abs2 + abs( TFview ) .^ 2;

    sum_conjT_exwv = sum_conjT_exwv + conj( TFview ) .* phi( :, :, :, ss );
    
%     for pp = 1 : sol.probe.scpm.N
%         sum_conjT_exwv( :, :, pp ) = sum_conjT_exwv( :, :, pp ) + conj( TFview ) .* phi( :, :, pp, ss ); end
    
    
    
    
    
end

P = sum_conjT_exwv ./ ( sum_abs_T_abs2 + 1e-7 );

% toc

5;










%{
    
tic

    [ ind ] = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );

    TFview = sol.sample.TF( : );
    TFview = reshape( TFview( ind ), [ sol.sz.sz, sol.spos.N ]);

    tmp0 = permute( repmat( conj( TFview ), [ 1, 1, 1, sol.probe.scpm.N ] ), [ 1, 2, 4, 3 ] );
    
    tmp1 = sum( tmp0 .* phi , 4 ) ./ sum( 1e-7 + abs( tmp0 ) .^ 2, 4 );
    
    toc
    
    
    figure; imagesc( abs( P( :, :, 2 ))); colorbar
    figure; imagesc( abs( tmp1( :, :, 2 ))); colorbar
    figure; imagesc( abs( sol.probe.P( :, :, 2))); colorbar
    
    norm( tmp1(:) - sol.probe.P( : ), 'fro' )
    
    5;
    
%     clear( 'TFview', 'probe' )


%}






























% P = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     P( :, :, pp ) = sum_conjT_exwv( :, :, pp ) ./ ( sum_abs_T_abs2 + 1e-7 ); 
% %     P( :, :, pp ) = 0.5 * P( :, :, pp ) + 0.5 * sum_conjT_exwv( :, :, pp ) ./ ( sum_abs_T_abs2 + 1e-7 ); 
%     
% end




