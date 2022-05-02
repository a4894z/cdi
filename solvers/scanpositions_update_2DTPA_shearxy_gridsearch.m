function [ opt_shear_y, opt_shear_x ] = scanpositions_update_2DTPA_shearxy_gridsearch( rs, T0, phi, measD, nDeq0, spos_opt )

    

shear_x = gpuArray( single( linspace( -0.5, 0.5, 5 )));
shear_y = gpuArray( single( linspace( -0.5, 0.5, 5 )));

L = gpuArray( zeros( length( shear_y ), length( shear_x ) , 'single' ) );

for sy = 1 : length( shear_y )
    
    for sx = 1 : length( shear_x )
        
        
        sysxT = gpuArray( [ [ 1, shear_y( sy ) ]; [ shear_x( sx ), 1 ] ]);
        
        
        rs_sysx = transpose( sysxT *  transpose ( rs ));

        ind = get_indices_2Dframes( rs_sysx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

%         I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        I_e =  sum( abs( psi ) .^ 2, 3 );
        
        L( sy, sx ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );


    end

end



[ ~, II ] = min( L(:) );

[ Ir, Ic ] = ind2sub( size( L ), II );

opt_shear_y = shear_y( Ir );
opt_shear_x = shear_y( Ic );


% figure; 
% imagesc( shear_x, shear_y, L ); 
% grid on
% 
% 
% 
% 5;

































%     ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     
% 
%     %===============
%     % search shear x
%     %===============
%     
%     rs_shearx = transpose( spos_opt.shear_x_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_shearx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_shearx = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_shearx = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
% 
%     grad_Lg_shearx = sign( Lg_0 - Lg_shearx );
% 
%     %========
%     
%     Lg_shearx_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_shear, 'single' ));
%     
%     for aa = 1 : spos_opt.Naalpha_shear
% 
%         s_x_alpha = spos_opt.shear_x + spos_opt.aalpha_shear( aa ) * grad_Lg_shearx;
%         
%         shear_x  = [ [ 1,          0 ]; ...
%                      [ s_x_alpha,  1 ] ];
%   
%         rs_shearx = transpose( shear_x * transpose( rs ));
% 
%         ind = get_indices_2Dframes( rs_shearx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%         I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         sqrt_I_e = sqrt( I_e );
% 
%         Lg_shearx_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;
% 
%     end
% 
%     [ ~, II ] = min( Lg_shearx_aalpha );
%     
%     shear_x = spos_opt.shear_x + spos_opt.aalpha_shear(  II ) * grad_Lg_shearx;
% 
% 
%     %===============
%     % search shear y
%     %===============
% 
%     ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     
%     %========
%     
%     rs_sheary = transpose( spos_opt.shear_y_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_sheary, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_sheary = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_shearx = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
% 
%     grad_Lg_sheary = sign( Lg_0 - Lg_sheary );
% 
%     %========
%     
%     Lg_sheary_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_shear, 'single' ));
%     
%     for aa = 1 : spos_opt.Naalpha_shear
% 
%         s_y_alpha = spos_opt.shear_y + spos_opt.aalpha_shear( aa ) * grad_Lg_sheary;
%         
%         shear_y  = [ [ 1,  s_y_alpha ]; ...
%                      [ 0,  1 ] ];
%   
%         rs_sheary = transpose( shear_y * transpose( rs ));
% 
%         ind = get_indices_2Dframes( rs_sheary, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%         I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         sqrt_I_e = sqrt( I_e );
% 
%         Lg_sheary_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;
% 
%     end
% 
%     [ ~, II ] = min( Lg_sheary_aalpha );
%     
%     shear_y = spos_opt.shear_y + spos_opt.aalpha_shear(  II ) * grad_Lg_sheary;


end

