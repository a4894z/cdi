function [ rs, rot ] = scanpositions_update_2DTPA_rotation( rs, T0, phi, measD, nDeq0, spos_opt )

% rs = fliplr( rs );

% nDeq0 = not( nDeq0 );
Nspos = size( rs, 1 );

rot_search = gpuArray( linspace( -30, 30, 61 ));

L = gpuArray( zeros( 1, length( rot_search ), 'single' ) );

for aa = 1 : length( rot_search )

    rotT = gpuArray( [ [ cosd( rot_search( aa ) ), -sind( rot_search( aa ) ) ]; ...
                       [ sind( rot_search( aa ) ), +cosd( rot_search( aa ) ) ] ] );     
        
                   
%     rotT = gpuArray( [ [ sind( rot_search( aa ) ), +cosd( rot_search( aa ) ) ]; ...
%                        [ cosd( rot_search( aa ) ), -sind( rot_search( aa ) ) ] ] );  
                   
    rs_rot = transpose( rotT * transpose( rs ));
            
    %========

    ind = get_indices_2Dframes( rs_rot, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = sum( abs( psi ) .^ 2, 3 );
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );

    L( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 )) / ( Nspos * spos_opt.rc );
%     L( aa ) = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( Nspos * spos_opt.rc );


end



[ ~, II ] = min( L );

rot = rot_search( II );

rotT = gpuArray( [ [ cosd( rot ), -sind( rot ) ] ; ...
                   [ sind( rot ), +cosd( rot ) ] ] );  
                   
rs = transpose( rotT * transpose( rs ));

% rs = fliplr( rs );


% rotT = gpuArray( [ [ sind( rot_search( II ) ), +cosd( rot_search( II ) ) ]; ...
%                    [ cosd( rot_search( II ) ), -sind( rot_search( II ) ) ] ] );     
% 
% rs_rot = transpose( rotT * transpose( rs ));



figure;
plot( rot_search, L )





%     ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( Nspos * spos_opt.rc );
% %     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( Nspos * spos_opt.rc );
%     
%     %======================
%     % search rotation angle
%     %======================
%     
%     
%     rs_scaley = transpose( spos_opt.scale_y_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_scaley = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( Nspos * spos_opt.rc );
% %     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( Nspos * spos_opt.rc );
% 
%     
    
    
    
    
    



end

