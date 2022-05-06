function [ rs, spos_opt ] = scanpositions_update_2DTPA_affine( rs, T0, phi, measD, nDeq0, spos_opt )

    %========
    
    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    I_e = sum( abs( psi ) .^ 2, 3 );

    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_0 = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos );     
        
    else
        
        L_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );  

    end

    %========

    rs_a = transpose( spos_opt.scale_y_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_a, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    I_e = sum( abs( psi ) .^ 2, 3 );

    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_a = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos );  
        
    else
        
        L_a = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );  

    end

    %========

    rs_b = transpose( spos_opt.shear_y_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_b, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    I_e = sum( abs( psi ) .^ 2, 3 );

    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_b = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos );     
        
    else
        
        L_b = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );  

    end
    
    %========

    rs_c = transpose( spos_opt.shear_x_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_c, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    I_e = sum( abs( psi ) .^ 2, 3 );

    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_c = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos );     
        
    else
        
        L_c = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );   

    end

    %========

    rs_d = transpose( spos_opt.scale_x_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_d, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    I_e = sum( abs( psi ) .^ 2, 3 );
    
    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_d = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos ); 
        
    else
        
        L_d = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );   

    end

    %========
    
    grad_L_a = L_0 - L_a;
    grad_L_b = L_0 - L_b;
    grad_L_c = L_0 - L_c;
    grad_L_d = L_0 - L_d;
    
    grad_L   = [ grad_L_a, grad_L_b, grad_L_c, grad_L_d ];
    grad_L = grad_L ./ ( 1e-7 + sqrt( sum( abs( grad_L ).^ 2 )));
    
    Laalpha_ALL = gpuArray( zeros( 1, spos_opt.Naalpha_affineT, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_affineT

        genT_aalpha = spos_opt.affineT + spos_opt.aalpha_affineT( aa ) * grad_L;
        
        genT_aalpha = [ genT_aalpha( 1 ), genT_aalpha( 2 ); genT_aalpha( 3 ), genT_aalpha( 4 ) ]; 
        
        rs_genT_aalpha = transpose( genT_aalpha * transpose( rs ));

        ind = get_indices_2Dframes( rs_genT_aalpha, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
%         I_e = squeeze( sum( abs( psi ) .^ 2, 3 ));
        I_e = sum( abs( psi ) .^ 2, 3 );
        
        if strcmp( spos_opt.noise_model, 'poisson' )

            Laalpha_ALL( aa ) = sum( sum( sum( nDeq0 .* ( I_e - measD .* log( I_e )) ))) / ( spos_opt.rc * spos_opt.Nspos );
     
        else
            
            Laalpha_ALL( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2 ))) / ( spos_opt.rc * spos_opt.Nspos );
 
        end

    end
    
    [ ~, II ] = min( Laalpha_ALL );
    spos_opt.affineT = spos_opt.affineT + spos_opt.aalpha_affineT( II ) * grad_L;
    
    genT_aalpha = [ spos_opt.affineT( 1 ), spos_opt.affineT( 2 ); spos_opt.affineT( 3 ), spos_opt.affineT( 4 ) ]; 
    
    rs = transpose( genT_aalpha * transpose( rs ));
 
end

