function [ rs ] = scanpositions_update_2DTPA( rs, T0, phi, measD, nDeq0, spos_opt )
    
    nDeq0 = not( nDeq0 );
    Nspos = size( rs, 1 );
    
    %========

    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e = sum( abs( psi ) .^ 2, 3 );
%     I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    
    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_0 = squeeze( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / spos_opt.rc;     
        
    else
        
        L_0 = squeeze( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 ))  / spos_opt.rc;

    end
    
    %========

    ind    = get_indices_2Dframes( rs + [ 1, 0 ], spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi_dy = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
    psi_dy = fft( fft( fftshift( fftshift( psi_dy, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e = sum( abs( psi_dy ) .^ 2, 3 );
%     I_e = squeeze( sum( abs( psi_dy ) .^ 2, 3 ) );

    if strcmp( spos_opt.noise_model, 'poisson' )
        
        L_dy = squeeze( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / spos_opt.rc;
        
    else
        
        L_dy = squeeze( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 ))  / spos_opt.rc;

    end
    
    %========

    ind    = get_indices_2Dframes( rs + [ 0, 1 ], spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi_dx = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
    psi_dx = fft( fft( fftshift( fftshift( psi_dx, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e = sum( abs( psi_dx ) .^ 2, 3 );
%     I_e = squeeze( sum( abs( psi_dx ) .^ 2, 3 ) );

    if strcmp( spos_opt.noise_model, 'poisson' )

        L_dx = squeeze( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / spos_opt.rc;

    else

        L_dx = squeeze( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 ))  / spos_opt.rc;

    end

    %========
    
    grad_L_y = L_0 - L_dy;
    grad_L_x = L_0 - L_dx;
    grad_L   = [ grad_L_y, grad_L_x ];
    grad_L_N = grad_L ./ sqrt( sum( abs( grad_L ).^ 2, 2 ));

    Laalpha     = gpuArray( zeros( Nspos, spos_opt.Naalpha_rs, 'single' ));
    Laalpha_ALL = gpuArray( zeros( 1, spos_opt.Naalpha_rs, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_rs

        ind = get_indices_2Dframes( rs + spos_opt.aalpha_rs( aa ) * grad_L_N, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

        I_e = sum( abs( psi ) .^ 2, 3 );
%         I_e = squeeze( sum( abs( psi ) .^ 2, 3 )) ;
        
        if strcmp( spos_opt.noise_model, 'poisson' )
            
            Laalpha( :, aa ) = squeeze( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )) , 1 ), 2 )) / spos_opt.rc;   
            Laalpha_ALL( aa ) = sum( Laalpha( :, aa ) ) / Nspos;

        else

            Laalpha( :, aa ) = squeeze( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 )) / spos_opt.rc;
            Laalpha_ALL( aa ) = sum( Laalpha( :, aa ) ) / Nspos; 
            
        end

    end

    %========
    
%     Laalpha = Laalpha - min( Laalpha, [], 2 );
%     Laalpha = Laalpha ./ max( Laalpha, [], 2 );

    %========

    if strcmp( spos_opt.optimize_rs_GD, 'individual' )
        
        [ ~, II ] = min( Laalpha, [], 2 );   
        rs = rs + transpose( spos_opt.aalpha_rs( II )) .* grad_L_N;
       
    else
        
        [ ~, II ] = min( Laalpha_ALL );
        rs = rs + spos_opt.aalpha_rs( II ) * grad_L_N;

    end

end

