function ptycho2DTPA_mkimg_meas_metric( sol, expt )

    
    tmp0 = sol.metrics.meas_gauss_intensity;
    tmp1 = sol.metrics.meas_gauss_magnitude;
    tmp2 = sol.metrics.meas_poiss;
    
    tmp0 = tmp0 - min( tmp0 );
    tmp1 = tmp1 - min( tmp1 );
    tmp2 = tmp2 - min( tmp2 );
    
    tmp0 = tmp0 / max( abs( tmp0 ));
    tmp1 = tmp1 / max( abs( tmp1 ));
    tmp2 = tmp2 / max( abs( tmp2 ));
    
    figure( 666 ); 
    set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     
    
    subplot(141)
    
    hold on
    plot( sol.it.mtot, log10( 1 + 10^3 * tmp0 ), '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )
    plot( sol.it.mtot, log10( 1 + 10^3 * tmp1 ), '-', 'Linewidth', 3, 'color', [ 1, 0, 0 ] )    
    plot( sol.it.mtot, log10( 1 + 10^3 * tmp2 ), '-', 'Linewidth', 3, 'color', [ 0, 1, 0 ] )
    hold off
    
    grid on
    
    title('log$_{10}( 1 + 10^3$ Normalized Metrics $ \in [ 0, 1 ])$', 'Interpreter', 'latex' )
    
    legend( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \vert I^m_{ \mathbf{r}, s } -  I^e_{ \mathbf{r}, s } \vert^2 $',                       ...
            '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \big\vert \sqrt{I^m_{ \mathbf{r}, s }} -  \sqrt{I^e_{ \mathbf{r}, s }} \big\vert^2 $', ...
            '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} I^e_{ \mathbf{r}, s } - I^m_{ \mathbf{r}, s } log( I^e_{ \mathbf{r}, s } ) $',         ...
            'FontSize', 14, ...
            'Interpreter','latex');

%     export_fig( 'meas_metric_gauss_poiss_normalized_zero_one.jpg', '-r90.0' )
%     close 666;

    %========================
    % measurement cost metric 
    %========================

%     figure( 666 ); 
%     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     
%     
    
    subplot(142)
    semilogy( sol.it.mtot, sol.metrics.meas_gauss_intensity, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )

    grid on
%     title( '$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} - \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 14,       ...
%            'Interpreter', 'latex' );

    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \vert I^m_{ \mathbf{r}, s } -  I^e_{ \mathbf{r}, s } \vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'Interpreter', 'latex' );
       
    
    
    subplot(143)
    semilogy( sol.it.mtot, sol.metrics.meas_gauss_magnitude, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )

    grid on
%     title( '$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} - \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F $', ...
%            'FontWeight', 'bold', ...
%            'FontSize', 14,       ...
%            'Interpreter', 'latex' );

    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} \big\vert \sqrt{I^m_{ \mathbf{r}, s }} -  \sqrt{I^e_{ \mathbf{r}, s }} \big\vert^2 $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'Interpreter', 'latex' );
       
%     export_fig( 'meas_metric_gauss.jpg', '-r90.0' )
%     close 666;
        
    
%     figure( 666 ); 
%     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     
    subplot(144)
    semilogy( sol.it.mtot, sol.metrics.meas_poiss, '-', 'Linewidth', 3, 'color', [ 0, 0, 0 ] )

    grid on
    title( '$ \frac{1}{N_s} \sum_s \sum_{\mathbf{r}} I^e_{ \mathbf{r}, s } - I^m_{ \mathbf{r}, s } log( I^e_{ \mathbf{r}, s } ) $', ...
           'FontWeight', 'bold', ...
           'FontSize', 14,       ...
           'Interpreter', 'latex' );

    export_fig( 'meas_metric_gauss_poiss.jpg', '-r90.0' )
    close 666;
        
    
end

