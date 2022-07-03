function [ Psi, poissonexwv, metrics ] = exitwave_update_2DTPA_poisson_v2( phi,       ...
                                                                          T0,         ...
                                                                          ind,        ...
                                                                          batch_indx, ...
                                                                          sz,         ...
                                                                          Nspos,      ...
                                                                          Nscpm,      ...
                                                                          poissonexwv, ...
                                                                          rc,         ... 
                                                                          sqrt_rc,    ...
                                                                          I_m,        ...
                                                                          I_m_eq0,    ...
                                                                          measLPF,    ...
                                                                          epoch,      ...
                                                                          collect_metrics )
                                                                     
    metrics = [];    

    % metrics.gauss = [];
    % metrics.poiss = [];

    %========================================
    % get sample frames at each scan position
    %========================================

    T = reshape( T0( ind ), [ sz, 1, Nspos ]);

    %===========================
    % form exitwaves using 2DTPA
    %===========================

    Psi = T .* phi;

    %=====================================
    % propagate exitwave modes to detector
    %=====================================

    Psi = fft( fft( fftshift( fftshift( Psi, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
    Psi = reshape( Psi, [ rc, Nscpm, Nspos ] );
    Psi = permute( Psi, [ 1, 3, 2 ] );

    %==================================================================
    % intensity for current sample/probe solutions at measurement plane
    %==================================================================

    abs2_Psi = abs( Psi ) .^ 2;

    I_e = sum( abs2_Psi, 3 ); % sum over scpm
    
    %========

    xi = 1 - I_m ./ I_e;

    % chknan = isnan( xi );
    % chkinf = isinf( xi );
    % 
    % if any( chknan( : ) ), xi( chknan ) = 0; end
    % if any( chkinf( : ) ), xi( chkinf ) = 0; end
    
    %========
    
%     grad_Psi = Psi .* xi;
    
    %=================================
    % Collect algo performance metrics
    %=================================
    
    if collect_metrics 

        meas_mask = not( I_m_eq0 ); % FIX: !!!! pass the already not version !!!!

        metrics.gauss_intensity = sum( sum( sum( meas_mask .* abs( I_m - I_e ) .^ 2                 )));
        metrics.gauss_magnitude = sum( sum( sum( meas_mask .* abs( sqrt( I_m ) - sqrt( I_e ) ) .^ 2 )));
        metrics.poiss           = sum( sum( sum( meas_mask .* ( I_e - I_m .* log( I_e ))            )));

        I_m_over_I_e = 1 - xi;
          
        metrics.grad_gauss_intensity = sum( sum( sum( meas_mask .* sum( abs( 2 * ( I_e - I_m )           .* Psi ) .^ 2, 3 ) )));
        metrics.grad_gauss_magnitude = sum( sum( sum( meas_mask .* sum( abs( ( 1 - sqrt( I_m_over_I_e )) .* Psi ) .^ 2, 3 ) )));
        metrics.grad_poiss           = sum( sum( sum( meas_mask .* sum( abs( ( 1 - I_m_over_I_e )        .* Psi ) .^ 2, 3 ) ))); 

    end

    %================================================
    % compute the step length by 1st order optimality
    %================================================
    
%     if ( strcmp( poissonexwv.steplength_update_type, 'use_exact_linesearch' ) && (( mod( epoch, poissonexwv.steplength_update_freq ) == 0 ) || ( epoch == 1 )))
% 
%         poissonexwv.alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi,                     ...
%                                                                                   abs2_Psi,               ...
%                                                                                   poissonexwv.alpha_test, ...
%                                                                                   I_e,                    ...
%                                                                                   I_m,                    ...
%                                                                                   Nspos,                  ...
%                                                                                   Nscpm );
%         
%         poissonexwv.alpha_prev( batch_indx, : ) = squeeze( poissonexwv.alpha_opt );
%     
%     elseif strcmp( poissonexwv.steplength_update_type, 'use_sign_test' ) && (( mod( epoch, poissonexwv.steplength_update_freq ) == 0 ) || ( epoch == 1 ))
% 
%         alpha_test = poissonexwv.alpha_prev( batch_indx, : );
%         alpha_test = reshape( alpha_test, [ 1, size( alpha_test ) ] );
%                                                                                                                                    
%         poissonexwv.alpha_opt = poisson_steplength_signtest_vs_minibatch( xi,         ...
%                                                                           abs2_Psi,   ...
%                                                                           alpha_test, ...
%                                                                           I_e,        ...
%                                                                           I_m,        ...
%                                                                           poissonexwv.delta_alpha_signtest );
% 
%         poissonexwv.alpha_prev( batch_indx, : ) = squeeze( poissonexwv.alpha_opt );
%         
%     elseif isfield( poissonexwv, 'alpha_prev' )
%         
%         poissonexwv.alpha_opt = reshape( poissonexwv.alpha_prev( batch_indx, : ), [ 1, Nspos, Nscpm ] );
%         
%     else
%         
%         error( '!!!!! Can''t Compute Optimal Step Length For Poisson Exitwave Update !!!!!' )
% 
%     end
%     
    
    
    
    
   if ( strcmp( poissonexwv.steplength_update_type, 'use_exact_linesearch' ) && (( mod( epoch, poissonexwv.steplength_update_freq ) == 0 ) || ( epoch == 1 )))

        alpha_prev = poissonexwv.alpha_prev( batch_indx, : );
        alpha_prev = reshape( alpha_prev, [ 1, size( alpha_prev ) ] );
        
        
        poissonexwv.alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi,                     ...
                                                                                  abs2_Psi,               ...
                                                                                  poissonexwv.alpha_test, ...
                                                                                  I_e,                    ...
                                                                                  I_m,                    ...
                                                                                  I_m_eq0,                ...
                                                                                  Nspos,                  ...
                                                                                  Nscpm,                  ...
                                                                                  alpha_prev );
        
        poissonexwv.alpha_prev( batch_indx, : ) = squeeze( poissonexwv.alpha_opt );
        
        %{
        
        figure; plot(squeeze( poissonexwv.alpha_opt )); 
        figure; plot(poissonexwv.alpha_prev( batch_indx, : ));
        figure; plot(poissonexwv.alpha_prev( batch_indx, : ) - squeeze( poissonexwv.alpha_opt ));
        
        %}
        
   else
       
        alpha_test = poissonexwv.alpha_prev( batch_indx, : );
        alpha_test = reshape( alpha_test, [ 1, size( alpha_test ) ] );
                                                                                                                                   
%         poissonexwv.alpha_opt = poisson_steplength_signtest_vs_minibatch( xi,         ...
%                                                                           abs2_Psi,   ...
%                                                                           alpha_test, ...
%                                                                           I_e,        ...
%                                                                           I_m,        ...
%                                                                           I_m_eq0,    ...
%                                                                           poissonexwv.delta_alpha_signtest );

                                               
                            
                                                                      
                                                                      

        poissonexwv.alpha_opt = poisson_steplength_recursive_vs_minibatch( xi,                          ...
                                                                           abs2_Psi,                    ...
                                                                           alpha_test,                  ...
                                                                           I_e,                         ...
                                                                           I_m,                         ...
                                                                           I_m_eq0,                     ...
                                                                           gpuArray( single( 0.5 )) );
                               


        poissonexwv.alpha_prev( batch_indx, : ) = squeeze( poissonexwv.alpha_opt );

   end
    
   
   
   
    % poisson_steplength_recursive_vs_minibatch

    
    
    
    
    %=======================================================================
    % apply the step length to make exitwaves satisfy measurement constraint
    %=======================================================================
    
    


    Psi = reshape( not( I_m_eq0 ) .* Psi .* ( 1 - poissonexwv.alpha_opt .* xi ) + Psi .* I_m_eq0 * 0.98, [ sz, Nspos, Nscpm ] );
%     Psi = reshape( Psi .* ( 1 - poissonexwv.alpha_opt .* xi ), [ sz, Nspos, Nscpm ] );
    
%     Psi = reshape( Psi - alpha_opt .* grad_Psi, [ sz, Nspos, Nscpm ] );





    


    %{
    
    Psi_missingdata = Psi .* I_m_eq0;
    Psi_missingdata_mat = reshape( Psi .* I_m_eq0, [ sz, Nspos, Nscpm ] );
    figure; imagesc( fftshift( log10( 1 + abs( Psi_missingdata_mat(:,:,22, 3) ) )))


    Psi_up = reshape( Psi .* ( 1 - poissonexwv.alpha_opt .* xi ), [ sz, Nspos, Nscpm ] );
    
    figure; imagesc( fftshift( log10( 1 + abs( Psi_up(:,:,22, 3) + Psi_missingdata_mat(:,:,22, 3) ) )))
    

    tmp0 = reshape( Psi, [ sz, Nspos, Nscpm ] );
    tmp1 = reshape( Psi.* ( 1 - poissonexwv.alpha_opt .* xi ), [ sz, Nspos, Nscpm ] );
    tmp2 = reshape( Psi .* ( 1 - poissonexwv.alpha_opt .* xi + 1.00 * I_m_eq0 ), [ sz, Nspos, Nscpm ] );
    tmp3 = reshape( not( I_m_eq0 ) .* Psi .* ( 1 - poissonexwv.alpha_opt .* xi ) + Psi .* I_m_eq0, [ sz, Nspos, Nscpm ] );


    figure; imagesc( fftshift( log10( 1 + abs( tmp0(:,:,22, 3) ) )))
    figure; imagesc( fftshift( log10( 1 + abs( tmp1(:,:,22, 3) ) )))
    figure; imagesc( fftshift( log10( 1 + abs( tmp2(:,:,22, 3) ) )))
    figure; imagesc( fftshift( log10( 1 + abs( tmp3(:,:,22, 3) ) )))





    %}





%     Psi = fftshift( fftshift( ifft( ifft( Psi .* measLPF, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    Psi = fftshift( fftshift( ifft( ifft( Psi, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;

    Psi = permute( Psi, [ 1, 2, 4, 3 ]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi, abs2_Psi, alpha_test, I_e, I_m, I_m_eq0, Nspos, Nscpm, alpha_prev  )

    [ lhs_steplength_eqn, rhs_steplength_eqn ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_test, I_e, I_m, I_m_eq0 );
 
    %========

    f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );

    [ ~, II ] = min( f_eq_0, [], 3 );

    alpha_opt = gpuArray.zeros( Nspos, Nscpm, 'single' );

    for pp = 1 : Nscpm

        alpha_opt( :, pp ) = alpha_test( 1, 1, pp, II( :, pp ) );

    end

    alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ] );


    
    
    
    
    
    
    
    
    
    
    
    %{
    
    sign_lhs_minus_rhs   = sign( lhs_steplength_eqn - rhs_steplength_eqn );
    delta_alpha_signtest = gpuArray( single( 0.05 ));
    
    tic
    alpha_STopt   = poisson_steplength_signtest_vs_minibatch( xi, abs2_Psi, alpha_prev, I_e, I_m, I_m_eq0, delta_alpha_signtest );
    toc
    
    tic
    alpha_RCopt_p = poisson_steplength_recursive_vs_minibatch( xi, abs2_Psi, alpha_prev, I_e, I_m, I_m_eq0, gpuArray( single( 0.5 )) );
    toc
    
    tic
    alpha_RCopt   = poisson_steplength_recursive_vs_minibatch_ptychoshelves( xi, alpha_prev, I_e, I_m, I_m_eq0, gpuArray( single( 0.5 )) );
    toc
    
    kk = 1;
    

    for pp = 1 : Nscpm

        %========

        figure( 667 );
        subplot( 1, double(gather( Nscpm )), kk )

        imagesc( squeeze( alpha_test( 1, 1, pp, : )), 1 : Nspos, squeeze( log10( 1 + f_eq_0( :, pp, : ) )))

        hold on
        plot( alpha_opt( 1, :, pp ), 1 : Nspos, 'rx', 'MarkerSize', 12, 'linewidth', 2 );
        plot( alpha_RCopt_p( 1, :, pp ), 1 : Nspos, '-o', 'MarkerSize', 6, 'color', [ 0.9, 0.0, 0.7 ], 'linewidth', 2 );
%         plot( alpha_STopt( 1, :, pp ), 1 : Nspos, '-o', 'MarkerSize', 6, 'color', [ 0.0, 0.7, 0.0 ], 'linewidth', 2 );
        plot( alpha_prev( 1, :, pp ), 1 : Nspos, 'bs', 'MarkerSize', 12, 'linewidth', 2 );
        plot( alpha_RCopt, 1 : Nspos, 'o', 'MarkerSize', 6, 'color', [ 0.9, 0.7, 0.0 ], 'linewidth', 2 );
        hold off

        grid on
        colormap bone
        set( gca, 'fontsize', 12 )
        xlabel('\alpha')
        ylabel('spos index')

        %========
        
        figure( 666	);
        subplot( 1, double(gather( Nscpm )), kk )

        imagesc( squeeze( alpha_test( 1, 1, pp, : )), 1 : Nspos, squeeze( sign_lhs_minus_rhs( :, pp, : ) ))
        
        hold on
        plot( alpha_opt( 1, :, pp ), 1 : Nspos, 'rx', 'MarkerSize', 12, 'linewidth', 2 );
        plot( alpha_RCopt_p( 1, :, pp ), 1 : Nspos, '-o', 'MarkerSize', 6, 'color', [ 0.9, 0.0, 0.7 ], 'linewidth', 2 );
%         plot( alpha_STopt( 1, :, pp ), 1 : Nspos, '-o', 'MarkerSize', 6, 'color', [ 0.0, 0.7, 0.0 ], 'linewidth', 2 );
        plot( alpha_prev( 1, :, pp ), 1 : Nspos, 'bs', 'MarkerSize', 12, 'linewidth', 2 );
        plot( alpha_RCopt, 1 : Nspos, 'o', 'MarkerSize', 6, 'color', [ 0.9, 0.7, 0.0 ], 'linewidth', 2 );
        hold off
        
        
        grid on
        colormap bone
        set( gca, 'fontsize', 12 )
        xlabel('\alpha')
        ylabel('spos index')

        
        
        

        %========
        
        kk = kk + 1;

        drawnow
        

    end

    %}
    
    
    
    
    
    
end

%====================================================================================================================================================

function alpha_opt = poisson_steplength_signtest_vs_minibatch( xi, abs2_Psi, alpha_prev, I_e, I_m, I_m_eq0, delta_alpha_signtest )

    [ lhs_steplength_eqn, rhs_steplength_eqn ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_prev, I_e, I_m, I_m_eq0 );
 
    sign_lhs_minus_rhs = sign( lhs_steplength_eqn - rhs_steplength_eqn );

    alpha_opt = squeeze( alpha_prev ) - delta_alpha_signtest * sign_lhs_minus_rhs;
    
    alpha_opt( alpha_opt < 0 ) = 0;
        
    %{
    
    Z1 = squeeze( alpha_test );
    Z2 = squeeze( alpha_opt );
    
    figure;
    plot( Z1(:,5), 'o')
    hold on
    plot( Z2(:,5), 'x')
    hold off

    %}
    
    alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ]);
    
end

%====================================================================================================================================================

function alpha_opt_p = poisson_steplength_recursive_vs_minibatch( xi, abs2_Psi, alpha_prev, I_e, I_m, I_m_eq0, w )

    I_m_eq0 = not( I_m_eq0 );
    
    %========
    
    xi_abs_Psi2        = xi .* abs2_Psi;
    xi_alpha_minus_one = xi .* alpha_prev - 1;
        
    numer = I_m .* xi_alpha_minus_one;
    denom = abs2_Psi .* xi_alpha_minus_one .^ 2 + I_e - abs2_Psi;
    
    %========
        
%     numer = sum( xi_abs_Psi2 .* ( 1 + numer ./ denom ), 1 );        % sum over image pixels
%     denom = sum( ( xi .^ 2 ) .* abs2_Psi, 1 );                      % sum over image pixels

    numer = sum( I_m_eq0 .* xi_abs_Psi2 .* ( 1 + numer ./ denom ), 1 );        % sum over image pixels
    denom = sum( I_m_eq0 .* ( xi .^ 2 ) .* abs2_Psi, 1 );                      % sum over image pixels
    
    %========
    
    alpha_opt_p = numer ./ denom;
    
    alpha_opt_p = w * alpha_opt_p + ( 1 - w ) * alpha_prev;
    
    alpha_opt_p( alpha_opt_p < 0 ) = 0;
    
end

%====================================================================================================================================================

function alpha_opt = poisson_steplength_recursive_vs_minibatch_ptychoshelves( xi, alpha_prev, I_e, I_m, I_m_eq0, w )

    I_m_eq0 = not( I_m_eq0 );
    
    %========
    
    alpha_opt = alpha_prev( 1, :, end );

    for ii = 1 : 2

        xi_alpha_minus_one = xi .* alpha_opt - 1;

        denom = ( xi .^ 2 ) .* I_e;

        numer = xi .* ( I_e + I_m ./ ( 1e-6 + xi_alpha_minus_one ));

    %     alpha_opt = squeeze( sum( numer, 1 ) ./  sum( denom, 1 ));

%         alpha_opt = w * squeeze( sum( numer, 1 ) ./  sum( denom, 1 )) + ( 1 - w ) * alpha_opt;
        alpha_opt = w * squeeze( sum( I_m_eq0 .* numer, 1 ) ./  sum( I_m_eq0 .* denom, 1 )) + ( 1 - w ) * alpha_opt;
        
    %     alpha_opt = max( min( alpha_opt, 1 ), 0 );

    %     z = ( alpha_opt > 1 );
    %     alpha_opt( z ) = alpha_prev( 1, z, end );

        z = ( alpha_opt < 0 );
        alpha_opt( z ) = alpha_prev( 1, z, end );

        z = isnan( alpha_opt );
        alpha_opt( z ) = alpha_prev( 1, z, end );

    
    end
    
    

end

%====================================================================================================================================================

function [ lhs, rhs ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_test, I_e, I_m, I_m_eq0 )
    
    I_m_eq0 = not( I_m_eq0 );
    
    %========
    
    xi_abs_Psi2        = xi .* abs2_Psi;
    xi_alpha_minus_one = xi .* alpha_test - 1;

    %========

%     lhs = squeeze( sum( xi_abs_Psi2 .* xi_alpha_minus_one, 1 ));    % sum over image pixels
    lhs = squeeze( sum( xi_abs_Psi2 .* xi_alpha_minus_one .* I_m_eq0, 1 ));    % sum over image pixels
    
    %========

    numer = I_m .* xi_alpha_minus_one;
    denom = abs2_Psi .* xi_alpha_minus_one .^ 2 + I_e - abs2_Psi;
    
%     rhs = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));    % sum over image pixels
    rhs = squeeze( sum( I_m_eq0 .* xi_abs_Psi2 .* ( numer ./ denom ), 1 ));    % sum over image pixels
    
end

%====================================================================================================================================================

