function [ sol ] = ptycho2DTPA_mkimg_meas_metric( sol, expt )

    figure( 666 ); 
    set( gcf, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )     

    if ~isfield( sol.metrics, 'poisson_offset' )
        
        I_m                        = expt.meas.D .^ 2;
        log_I_m                    = log( I_m );
        log_I_m( isinf( log_I_m )) = 0;

        sol.metrics.poisson_offset = ( I_m - I_m .* log_I_m );
        sol.metrics.poisson_offset = sum( sol.metrics.poisson_offset(:) ) / size( I_m, 3 );

        clear( 'I_m', 'log_I_m' )

    end

    tmp0 = sol.metrics.meas_gauss_intensity;
    tmp1 = sol.metrics.meas_gauss_magnitude;
    tmp2 = sol.metrics.meas_poiss - sol.metrics.poisson_offset;

    tmp3 = sol.metrics.grad_meas_gauss_intensity;
    tmp4 = sol.metrics.grad_meas_gauss_magnitude;
    tmp5 = sol.metrics.grad_meas_poiss;
    
    tmp0 = tmp0 / max( abs( tmp0 ));
    tmp1 = tmp1 / max( abs( tmp1 ));
    tmp2 = tmp2 / max( abs( tmp2 ));

    tmp3 = tmp3 / max( abs( tmp3 ));
    tmp4 = tmp4 / max( abs( tmp4 ));
    tmp5 = tmp5 / max( abs( tmp5 ));
   
    %===========================================
    % GAUSSIAN INTENSITY measurement cost metric 
    %===========================================

    subplot(231)
    semilogy( sol.it.mtot, tmp0, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Gaussian Intensity Metric')
    
    ylim( sol.metrics.metrics_ylim )
    
    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )

    legend( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \vert I^m_{ \mathbf{q}, s } -  I^e_{ \mathbf{q}, s } \vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
    %===========================================
    % GAUSSIAN MAGNITUDE measurement cost metric 
    %===========================================
    
    subplot(232)
    semilogy( sol.it.mtot, tmp1, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Gaussian Magnitude Metric')
    
    ylim( sol.metrics.metrics_ylim )
    
    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
    
    legend( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \big\vert \sqrt{I^m_{ \mathbf{q}, s }} -  \sqrt{I^e_{ \mathbf{q}, s }} \big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
    %================================
    % POISSON measurement cost metric 
    %================================

    subplot(233)
    semilogy( sol.it.mtot, tmp2, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Poisson Metric')
    
    ylim( sol.metrics.metrics_ylim )
    
    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
    
    legend( '$ -C_0 + \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} I^e_{ \mathbf{q}, s } - I^m_{ \mathbf{q}, s } log( I^e_{ \mathbf{q}, s } ) $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
    title( '$ C_0 = \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} I^m_{ \mathbf{q}, s } - I^m_{ \mathbf{q}, s } log( I^m_{ \mathbf{q}, s } ) $', ...
            'FontWeight', 'bold', ...
            'FontSize', 14,       ...
            'Interpreter', 'latex' );
    
    %===============================================================
    % Norm of Gradient of GAUSSIAN INTENSITY measurement cost metric 
    %===============================================================
    
    subplot(234)
    semilogy( sol.it.mtot, tmp3, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Gaussian Intensity, First Order Optimality')
    
    ylim( sol.metrics.grad_metrics_ylim )
    
    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )

    legend( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \big\vert  2 \tilde{\Psi}_{\mathbf{q},s,p} \big( I^e_{ \mathbf{q}, s } -  I^m_{ \mathbf{q}, s } \big) \big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
    %===============================================================
    % Norm of Gradient of GAUSSIAN MAGNITUDE measurement cost metric 
    %===============================================================

    subplot(235)
    semilogy( sol.it.mtot, tmp4, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Gaussian Magnitude, First Order Optimality')
    
    ylim( sol.metrics.grad_metrics_ylim )

    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )

    legend( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \Big\vert \tilde{\Psi}_{\mathbf{q},s,p} \Big( 1 - \sqrt{\frac{I^m_{ \mathbf{q}, s } }{I^e_{ \mathbf{q}, s } }} \Big) \Big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
    %====================================================
    % Norm of Gradient of POISSON measurement cost metric 
    %====================================================

    subplot(236)
    semilogy( sol.it.mtot, tmp5, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    xlabel('Epoch')
    ylabel('Poisson, First Order Optimality')
    
    ylim( sol.metrics.grad_metrics_ylim )

    grid on
    set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )

    legend( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \Big\vert \tilde{\Psi}_{\mathbf{q},s,p} \Big( 1 - \frac{I^m_{ \mathbf{q}, s } }{I^e_{ \mathbf{q}, s } } \Big) \Big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'location', sol.metrics.legend_loc, ...
           'Interpreter', 'latex' );
       
       
       
       
%     
%     subplot(141)
%     
%     hold on
%     plot( sol.it.mtot, log10( 1 + 10^3 * tmp0 ), '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
%     plot( sol.it.mtot, log10( 1 + 10^3 * tmp1 ), '-', 'Linewidth', 3, 'color', [ 1, 0, 0 ] )    
%     plot( sol.it.mtot, log10( 1 + 10^3 * tmp2 ), '-', 'Linewidth', 3, 'color', [ 0, 1, 0 ] )
%     hold off
%     
%     grid on
%     
%     title('log$_{10}( 1 + 10^3$ Normalized Metrics $ \in [ 0, 1 ])$', 'Interpreter', 'latex' )
%     
%     legend( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \vert I^m_{ \mathbf{r}, s } -  I^e_{ \mathbf{r}, s } \vert^2 $',                       ...
%             '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \big\vert \sqrt{I^m_{ \mathbf{r}, s }} -  \sqrt{I^e_{ \mathbf{r}, s }} \big\vert^2 $', ...
%             '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} I^e_{ \mathbf{r}, s } - I^m_{ \mathbf{r}, s } log( I^e_{ \mathbf{r}, s } ) $',         ...
%             'FontSize', 14, ...
%             'Interpreter','latex');
    
%     %===========================================
%     % GAUSSIAN INTENSITY measurement cost metric 
%     %===========================================
% 
%     subplot(231)
%     semilogy( sol.it.mtot, sol.metrics.meas_gauss_intensity, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
% 
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
% 
%     title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \vert I^m_{ \mathbf{q}, s } -  I^e_{ \mathbf{q}, s } \vert^2 $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
%        
%     %===========================================
%     % GAUSSIAN MAGNITUDE measurement cost metric 
%     %===========================================
%     
%     subplot(232)
%     semilogy( sol.it.mtot, sol.metrics.meas_gauss_magnitude, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
% 
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
%     
%     title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} \big\vert \sqrt{I^m_{ \mathbf{q}, s }} -  \sqrt{I^e_{ \mathbf{q}, s }} \big\vert^2 $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
%        
%     %================================
%     % POISSON measurement cost metric 
%     %================================
% 
%     subplot(233)
%     semilogy( sol.it.mtot, sol.metrics.meas_poiss - sol.metrics.poisson_offset, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
%     
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
%     
%     title( '$ C_0 + \frac{1}{N_s} \sum_s \sum_{\mathbf{q}} I^e_{ \mathbf{q}, s } - I^m_{ \mathbf{q}, s } log( I^e_{ \mathbf{q}, s } ) $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
% 
%        
%     %===============================================================
%     % Norm of Gradient of GAUSSIAN INTENSITY measurement cost metric 
%     %===============================================================
% 
%     subplot(234)
%     semilogy( sol.it.mtot, sol.metrics.grad_meas_gauss_intensity, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
% 
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
% 
%     title( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \big\vert  2 \tilde{\Psi}_{\mathbf{q},s,p} \big( I^e_{ \mathbf{q}, s } -  I^m_{ \mathbf{q}, s } \big) \big\vert^2 $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
%        
%     %===============================================================
%     % Norm of Gradient of GAUSSIAN MAGNITUDE measurement cost metric 
%     %===============================================================
% 
%     subplot(235)
%     semilogy( sol.it.mtot, sol.metrics.grad_meas_gauss_magnitude, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
% 
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
% 
%     title( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \Big\vert \tilde{\Psi}_{\mathbf{q},s,p} \Big( 1 - \sqrt{\frac{I^m_{ \mathbf{q}, s } }{I^e_{ \mathbf{q}, s } }} \Big) \Big\vert^2 $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
%        
%     %====================================================
%     % Norm of Gradient of POISSON measurement cost metric 
%     %====================================================
% 
%     subplot(236)
%     semilogy( sol.it.mtot, sol.metrics.grad_meas_poiss, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
% 
%     grid on
%     set( gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'times' )
% 
%     title( '$ \frac{1}{N_s N_p} \sum_s \sum_{\mathbf{q}} \sum_p \Big\vert \tilde{\Psi}_{\mathbf{q},s,p} \Big( 1 - \frac{I^m_{ \mathbf{q}, s } }{I^e_{ \mathbf{q}, s } } \Big) \Big\vert^2 $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 16,       ...
%            'Interpreter', 'latex' );
%        
       
    %=========================================
    % Finalize and save figure of cost metrics
    %=========================================
    
    export_fig( 'meas_metric_gauss_poiss.jpg', '-r120.0' )
    close 666;
        
    
end

