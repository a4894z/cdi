function [ opt_scale_y, opt_scale_x ] = scanpositions_update_2DTPA_scalexy_gridsearch( sol, expt )




Nrr = 50;

% scalex_search = gpuArray( single( 1 + 1 * linspace( -0.10, 0.10, 3 ) ));
% scaley_search = gpuArray( single( 1 + 1 * linspace( -0.10, 0.10, 3 ) ));
scalex_search = gpuArray( single( 1 + 1 * linspace( -0.05, 0.05, 3 ) ));
scaley_search = gpuArray( single( 1 + 1 * linspace( -0.05, 0.05, 3 ) ));
% scalex_search = gpuArray( single( 1 + 1 * linspace( -0.025, 0.025, 3 ) ));
% scaley_search = gpuArray( single( 1 + 1 * linspace( -0.025, 0.025, 3 ) ));
% scalex_search = gpuArray( single( 1 + 1 * linspace( -0.01, 0.01, 3 ) ));
% scaley_search = gpuArray( single( 1 + 1 * linspace( -0.01, 0.01, 3 ) ));
% scalex_search = gpuArray( single( 1 + 1 * linspace( -0.005, 0.005, 3 ) ));
% scaley_search = gpuArray( single( 1 + 1 * linspace( -0.005, 0.005, 3 ) ));

gauss_magnitude_scalexy = gpuArray.zeros( length( scaley_search ), length( scalex_search ), Nrr, 'single' );
            


GPU.batch_indx = gpuArray( single( 1 : 5 : sol.spos.N ));
% GPU.batch_indx = gpuArray( single( randperm( sol.spos.N, round( 0.15 * sol.spos.N ) )));

GPU.Nspos = gpuArray( length( GPU.batch_indx ) );

GPU.batch_rs = gpuArray( sol.spos.rs( GPU.batch_indx, : ) );

GPU.meas     = gpuArray( reshape( expt.meas.D( :, :, GPU.batch_indx ), [ sol.GPU.sz, 1, GPU.Nspos ] ));
GPU.meas_eq0 = ( GPU.meas == 0 );

    
for scx = 1 : length( scalex_search )
    
    for scy = 1 : length( scaley_search )
                       
    dt = tic;
    
    GPU.TFvec = sol.GPU.TFvec;

    
    
    
    GPU.samsz = sol.GPU.samsz;
    GPU.samrc = sol.GPU.samrc;
    GPU.vs_r = sol.GPU.vs_r;
    GPU.vs_c = sol.GPU.vs_c;


% GPU.TFvec = padarray( reshape( GPU.TFvec, sol.GPU.samsz ), [ 50, 50 ], 1 );
% 
% GPU.samsz = gpuArray( single( size( GPU.TFvec )));
% GPU.samrc = GPU.samsz( 1 ) * GPU.samsz( 2 );
% GPU.vs_r = gpuArray( single( round( ( 0.5 * ( GPU.samsz(1) - sol.GPU.sz(1) ) + 1 ) : ( 0.5 * ( GPU.samsz(1) + sol.GPU.sz(1) )))));
% GPU.vs_c = gpuArray( single( round( ( 0.5 * ( GPU.samsz(2) - sol.GPU.sz(2) ) + 1 ) : ( 0.5 * ( GPU.samsz(2) + sol.GPU.sz(2) )))));
% 
% GPU.TFvec = GPU.TFvec( : );                 
%         




        GPU.rs_search( :, 2 ) = scalex_search( scx ) * GPU.batch_rs( :, 2 );
        GPU.rs_search( :, 1 ) = scaley_search( scy ) * GPU.batch_rs( :, 1 );

        %========
                        
        GPU.ind_new = gpuArray( uint32( get_indices_2Dframes( GPU.rs_search, GPU.samsz, GPU.vs_r, GPU.vs_c ))); 

        GPU.ind_new_offset = GPU.ind_new + gpuArray( uint32( GPU.samrc * ( 0 : 1 : ( GPU.Nspos - 1 ) )));

        
                      
    for rr = 1 : Nrr

%         GPU.rs_search( :, 2 ) = scalex_search( scx ) * GPU.batch_rs( :, 2 );
%         GPU.rs_search( :, 1 ) = scaley_search( scy ) * GPU.batch_rs( :, 1 );
% 
%         %========
%                         
%         GPU.ind_new = gpuArray( uint32( get_indices_2Dframes( GPU.rs_search, GPU.samsz, GPU.vs_r, GPU.vs_c ))); 
% 
%         GPU.ind_new_offset = GPU.ind_new + gpuArray( uint32( GPU.samrc * ( 0 : 1 : ( GPU.Nspos - 1 ) )));

        %========

        % get exitwaves that satisfy the measurement constraints
        [ GPU.psi, metrics_rot ] = exitwave_update_2DTPA_gaussian( sol.GPU.phi,      ...            
                                                                   GPU.TFvec,    ...
                                                                   GPU.ind_new,      ...
                                                                   sol.GPU.sz,       ...
                                                                   GPU.Nspos,    ...
                                                                   sol.GPU.sqrt_rc,  ...
                                                                   GPU.meas,     ...
                                                                   GPU.meas_eq0, ...
                                                                   sol.GPU.measLPF,  ...
                                                                   logical( 1 ) );

        %========
        
        % get new sample    
        [ GPU.TFvec ] = rPIEupdate_batch_2DTPA_sample( GPU.psi,        ...
                                                       GPU.TFvec,      ...
                                                       sol.GPU.phi,        ...
                                                       GPU.ind_new_offset, ...
                                                       sol.GPU.rc,         ...
                                                       GPU.Nspos,      ...
                                                       sol.GPU.rPIE_alpha_T );  

        GPU.TFvec = modulus_limits_project( GPU.TFvec, sol.GPU.abs_TF_lim );
        
        %========
        
        gauss_magnitude_scalexy( scy, scx, rr ) = metrics_rot.gauss_magnitude / sol.GPU.Nspos;
%                         gauss_magnitude_shearxy( shy, shx, rr ) = metrics_rot.gauss_magnitude / sol.GPU.Nspos;
                         
    end
                    
    t = toc( dt );

    fprintf( [ num2str( [ scalex_search( scx ), scaley_search( scy ), t ], '( scx, scy ) = ( %0.3f, %0.3f ), t = %0.3f' ), '\n' ] )


    end
end
            

min_cost = min( gauss_magnitude_scalexy( : ));
max_cost = max( gauss_magnitude_scalexy( : ));

sz = size( gauss_magnitude_scalexy );

            
for rr = [ 1, 50 : 50 : Nrr ]
% for rr = 1 : 1 : Nrr
% for rr = round( linspace( 1, Nrr, 20 ))
                
    A          = gauss_magnitude_scalexy( :, :, rr );
    [ ~, II ]  = min( A(:) );
    [ ir, ic ] = ind2sub( [ sz( 1 ), sz( 2 ) ], II );
              
%                 figure;  
    set( figure, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
%                 imagesc( scalex_search, scaley_search, gauss_magnitude_scalexy( :, :, lose arr ), [ min_cost, max_cost ] )
    imagesc( scalex_search, scaley_search, gauss_magnitude_scalexy( :, :, rr ))  
    hold on
    plot( scalex_search( ic ), scaley_search( ir ), 'x', 'markersize', 30, 'linewidth', 2, 'color', [ 0.8, 0.0, 0.0 ] )
    hold off
    daspect([1 1 1]); 
    colorbar; 
    colormap turbo;                set( gca, 'ydir', 'normal' )
    title(num2str(rr, 'Epoch = %d'))
    grid on

%     export_fig( num2str( rr, 'scalexy_cost_%d.jpg' ), '-r120.0' )
%     rr
        
end
            
close all;


opt_scale_y = gather( scaley_search( ir ));
opt_scale_x = gather( scalex_search( ic ));

% scalex_search           = gather( scalex_search );
% scaley_search           = gather( scaley_search );
% gauss_magnitude_scalexy = gather( gauss_magnitude_scalexy );
            


5;



























% function [ opt_scale_y, opt_scale_x ] = scanpositions_update_2DTPA_scalexy_gridsearch( rs, T0, phi, measD, nDeq0, spos_opt )
% 
% scale_x = linspace( 0.8, 1.2, 5 );
% scale_y = linspace( 0.8, 1.2, 5 );
% 
% L = gpuArray( zeros( length( scale_y ), length( scale_x ) , 'single' ) );
% 
% for sy = 1 : length( scale_y )
%     
%     for sx = 1 : length( scale_x )
%         
%         
%         sysxT = gpuArray( [ [ scale_y( sy ), 0 ]; [ 0, scale_x( sx ) ] ]);
%         
%         
%         rs_sysx = transpose( sysxT *  transpose ( rs ));
% 
%         ind = get_indices_2Dframes( rs_sysx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
% %         I_e = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         I_e = sum( abs( psi ) .^ 2, 3 );
% 
%         L( sy, sx ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt( I_e ) ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
% 
%     end
% 
% end
% 
% 
% 
% [ ~, II ] = min( L(:) );
% 
% [ Ir, Ic ] = ind2sub( size( L ), II );
% 
% opt_scale_y = scale_y( Ir );
% opt_scale_x = scale_y( Ic );
% 
% % 
% % figure; 
% % imagesc( scale_x, scale_y, L )
% % grid on
% % 
% % 
% % 5;














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

