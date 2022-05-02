function [ opt_scale_y, opt_scale_x ] = scanpositions_update_2DTPA_scalexy_gridsearch( rs, T0, phi, measD, nDeq0, spos_opt )



scale_x = linspace( 0.8, 1.2, 5 );
scale_y = linspace( 0.8, 1.2, 5 );

L = gpuArray( zeros( length( scale_y ), length( scale_x ) , 'single' ) );

for sy = 1 : length( scale_y )
    
    for sx = 1 : length( scale_x )
        
        
        sysxT = gpuArray( [ [ scale_y( sy ), 0 ]; [ 0, scale_x( sx ) ] ]);
        
        
        rs_sysx = transpose( sysxT *  transpose ( rs ));

        ind = get_indices_2Dframes( rs_sysx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

%         I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        I_e = sum( abs( psi ) .^ 2, 3 );

        L( sy, sx ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );


    end

end



[ ~, II ] = min( L(:) );

[ Ir, Ic ] = ind2sub( size( L ), II );

opt_scale_y = scale_y( Ir );
opt_scale_x = scale_y( Ic );

% 
% figure; 
% imagesc( scale_x, scale_y, L )
% grid on
% 
% 
% 5;














%     %===============
%     % search scale y
%     %===============
%     
%     rs_scaley = transpose( spos_opt.scale_y_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_scaley = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
% 
%     grad_Lg_scaley = sign( Lg_0 - Lg_scaley );
% 
%     %========
%     
%     Lg_scaley_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_scale, 'single' ));
%     
%     for aa = 1 : spos_opt.Naalpha_scale
% 
%         scale_y_alpha = spos_opt.scale_y + spos_opt.aalpha_scale( aa ) * grad_Lg_scaley;
%         
%         rs_scaley = [ scale_y_alpha * rs( :, 1 ), rs( :, 2 ) ];
% 
%         ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%         I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         sqrt_I_e = sqrt( I_e );
% 
%         Lg_scaley_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;
% 
%     end
% 
%     [ ~, II ] = min( Lg_scaley_aalpha );
%     
%     scale_y = spos_opt.scale_y + spos_opt.aalpha_scale( II ) * grad_Lg_scaley;
% 
%     %===============
%     % search scale x
%     %===============
%     
%     rs_scalex = transpose( spos_opt.scale_x_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_scalex, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_scalex = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
% 
%     grad_Lg_scalex = sign( Lg_0 - Lg_scalex );
% 
%     %========
%     
%     Lg_scalex_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_scale, 'single' ));
%     
%     for aa = 1 : spos_opt.Naalpha_scale
% 
%         scale_x_alpha = spos_opt.scale_x + spos_opt.aalpha_scale( aa ) * grad_Lg_scalex;
%         
%         rs_scalexy = [ rs( :, 1 ), scale_x_alpha * rs( :, 2 ) ];
% 
%         ind = get_indices_2Dframes( rs_scalexy, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%         I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         sqrt_I_e = sqrt( I_e );
% 
%         Lg_scalex_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;
% 
%     end
% 
%     [ ~, II ] = min( Lg_scalex_aalpha );
%     
%     scale_x = spos_opt.scale_x + spos_opt.aalpha_scale( II ) * grad_Lg_scalex;    
    
end

