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
    
    if strcmp( poissonexwv.steplength_update_type, 'use_exact_linesearch' ) && (( mod( epoch, poissonexwv.steplength_update_freq ) == 0 ) || ( epoch == 1 ))

        poissonexwv.alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi,                     ...
                                                                                  abs2_Psi,               ...
                                                                                  poissonexwv.alpha_test, ...
                                                                                  I_e,                    ...
                                                                                  I_m,                    ...
                                                                                  Nspos,                  ...
                                                                                  Nscpm );
        
        poissonexwv.alpha_prev( batch_indx, : ) = squeeze( poissonexwv.alpha_opt );
    
        
        
        
        
        
        
        
        alpha_test = poissonexwv.alpha_prev( batch_indx, : );
        alpha_test = reshape( alpha_test, [ 1, size( alpha_test ) ] );
        
%         poissonexwv.alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi,                     ...
%                                                                                   abs2_Psi,               ...
%                                                                                   alpha_test, ...
%                                                                                   I_e,                    ...
%                                                                                   I_m,                    ...
%                                                                                   Nspos,                  ...
%                                                                                   Nscpm );
                                                                              
                                                                              
                                                                              
        alpha_opt = poisson_steplength_signtest_vs_minibatch( xi, abs2_Psi, alpha_test, I_e, I_m, Nspos, Nscpm );
        
        
    elseif strcmp( poissonexwv.steplength_update_type, 'use_sign_test' ) && (( mod( epoch, poissonexwv.steplength_update_freq ) == 0 ) || ( epoch == 1 ))

        alpha_test = poissonexwv.alpha_prev( batch_indx, : );
        alpha_test = reshape( alpha_test, [ 1, size( alpha_test ) ] );
        
        
    elseif isfield( poissonexwv, 'alpha_prev' )
        
        poissonexwv.alpha_opt = reshape( poissonexwv.alpha_prev( batch_indx, : ), [ 1, Nspos, Nscpm ] );
        
    else
        
        error( '!!!!! Can''t Compute Optimal Step Length For Poisson Exitwave Update !!!!!' )

    end

    %=======================================================================
    % apply the step length to make exitwaves satisfy measurement constraint
    %=======================================================================
    
%     tic
    Psi = reshape( Psi .* ( 1 - poissonexwv.alpha_opt .* xi ), [ sz, Nspos, Nscpm ] );
%     Psi = reshape( Psi - alpha_opt .* grad_Psi, [ sz, Nspos, Nscpm ] );
%     toc
    
%     Psi = fftshift( fftshift( ifft( ifft( Psi .* measLPF, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    Psi = fftshift( fftshift( ifft( ifft( Psi, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;

    Psi = permute( Psi, [ 1, 2, 4, 3 ]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha_opt = poisson_steplength_exact_linesearch_vs_minibatch( xi, abs2_Psi, alpha_test, I_e, I_m, Nspos, Nscpm )

    [ lhs_steplength_eqn, rhs_steplength_eqn ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_test, I_e, I_m );
 
    %========

    f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );

    [ ~, II ] = min( f_eq_0, [], 3 );

    alpha_opt = gpuArray.zeros( Nspos, Nscpm, 'single' );

    for pp = 1 : Nscpm

        alpha_opt( :, pp ) = alpha_test( 1, 1, pp, II( :, pp ) );

    end

    alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ] );


    
    
    
    %%{
    
    sign_lhs_minus_rhs = sign( lhs_steplength_eqn - rhs_steplength_eqn );
        
    kk = 1;

    for pp = 1 : Nscpm

        %========

        figure( 667 );
        subplot( 1, double(gather( Nscpm )), kk )

        imagesc( squeeze( alpha_test( 1, 1, pp, : )), 1 : Nspos, squeeze( log10( 1 + f_eq_0( :, pp, : ) )))

        hold on
        plot( alpha_opt( 1, :, pp ), 1 : Nspos, 'r.', 'MarkerSize', 5 );
        hold off

        grid on
        colormap bone
        xlabel('\alpha')
        ylabel('spos index')

        %========
        
        figure( 666	);
        subplot( 1, double(gather( Nscpm )), kk )
        

        
        imagesc( squeeze( alpha_test( 1, 1, pp, : )), 1 : Nspos, squeeze( sign_lhs_minus_rhs( :, pp, : ) ))
        
        hold on
        plot( alpha_opt( 1, :, pp ), 1 : Nspos, 'r.', 'MarkerSize', 5 );
        hold off
        
        grid on
        colormap bone
        xlabel('\alpha')
        ylabel('spos index')

        
        
        

        %========
        
        kk = kk + 1;

        drawnow
        
        
        

    end

    %}
    
    5;
    
end

%====================================================================================================================================================

function alpha_opt = poisson_steplength_signtest_vs_minibatch( xi, abs2_Psi, alpha_test, I_e, I_m, delta_alpha_signtest )

    [ lhs_steplength_eqn, rhs_steplength_eqn ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_test, I_e, I_m );
 
    sign_lhs_minus_rhs = sign( lhs_steplength_eqn - rhs_steplength_eqn );

    alpha_opt = squeeze( alpha_test ) - delta_alpha_signtest * sign_lhs_minus_rhs;
    
%     Z1 = squeeze( alpha_test );
%     Z2 = squeeze( alpha_opt );
%     
%     figure;
%     plot( Z1(:,5), 'o')
%     hold on
%     plot( Z2(:,5), 'x')
%     hold off

    alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ]);
    
end

%====================================================================================================================================================

function [ lhs, rhs ] = compute_lhs_rhs_optimal_poisson_steplength( xi, abs2_Psi, alpha_test, I_e, I_m )

    xi_abs_Psi2        = xi .* abs2_Psi;
    xi_alpha_minus_one = xi .* alpha_test - 1;

    %========

    lhs = squeeze( sum( xi_abs_Psi2 .* xi_alpha_minus_one, 1 ));    % sum over image pixels

    %========

    numer = I_m .* xi_alpha_minus_one;
    denom = abs2_Psi .* xi_alpha_minus_one .^ 2 + I_e - abs2_Psi;

    rhs = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));    % sum over image pixels

end
