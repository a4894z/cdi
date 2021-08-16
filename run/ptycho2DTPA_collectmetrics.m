function [ sol ] = ptycho2DTPA_collectmetrics( sol, expt )

    %=================================
    % get ready for array broadcasting
    %=================================

    meas_D    = reshape( expt.meas.D, [ expt.sz.sz, 1, expt.spos.N ] );
    meas_Deq0 = not( meas_D == 0 );

    %=========================================
    % compute exitwaves for all scan positions
    %=========================================

    spos.batch_indxsubset = 1 : sol.spos.N;
    spos.batch_rs         = sol.spos.rs( spos.batch_indxsubset, : );
    spos.frameindx        = get_indices_2Dframes( spos.batch_rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );

    %========

    TFv  = sol.sample.T( : );
    TF   = reshape( TFv( spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N ]);
    tmp0 = TF .* sol.probe.phi;
    tmp1 = fft( fft( fftshift( fftshift( tmp0, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;

    %===========================================================================
    % standard Gaussian noise metric with constant (ignored) stdev at all pixels
    %===========================================================================

    meas_residual = meas_Deq0 .* sqrt( sum( abs( tmp1 ) .^ 2, 3 )) - meas_D;
%     meas_residual = squeeze( sqrt( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 )));
    meas_residual = squeeze( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 ));
    
%     clear( 'tmp1', 'TF' )

%     %============================
%     % RAAR exitwave change metric
%     %============================
% 
%     psi = RAAR_GPU_arrays_hadamard_v2(  tmp0,            ...
%                                         sol.probe.phi,     ...
%                                         TFv,             ...
%                                         spos.frameindx,  ...
%                                         sol.sz.sz,       ...
%                                         sol.spos.N,      ...
%                                         sol.sz.sqrt_rc,  ...
%                                         meas_D,          ...
%                                         meas_Deq0,       ...
%                                         sol.measLPF,     ...
%                                         sol.RAAR.beta );
% 
% 
%     raar_exwv_change = abs( psi - tmp0 ) .^ 2;
% 
%     clear( 'tmp0', 'TFv', 'meas_D', 'meas_Deq0', 'spos' )

    %================
    % collect metrics
    %================

%     raar_exwv_change_batch = raar_exwv_change( sol.spos.batch_indxsubset );
%     raar_exwv_change_batch = raar_exwv_change_batch( : );

%     sol.metrics.exwv_change( sol.it.metr ) = sum( raar_exwv_change_batch ) / length( raar_exwv_change_batch );
    sol.metrics.meas_all( sol.it.metr ) = sum( meas_residual ) / length( meas_residual );
    
    %========

    % fprintf( [ '\n\n', num2str( [ kk, Nit, sol.it.epoch, sol.metrics.meas( sol.it.metr ), sol.metrics.exwv_change( sol.it.metr )], ...
    %             'iteration = %d / %d, iter total = %d, meas metric = %.2f, exit wave sample probe difference = %.2f' ), '\n\n' ]);

    sol.it.mtot( sol.it.metr ) = sol.it.epoch;

    %========

    figure( 666 ); 
    set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     

%     subplot( 2, 1, 1 ); 

    hold on

    semilogy( sol.it.mtot, sol.metrics.meas_all, '--', 'Linewidth', 4, 'color', [ 0, 0, 0 ] )

    % semilogy( sol.it.mtot, sol.metrics.meas, '-o', 'Linewidth', 2, 'Color', [0.8, 0, 0 ] ); 
    % semilogy( sol.it.mtot, sol.metrics.meas_IN, '-o', 'Linewidth', 2, 'Color', [0.0, 0.8, 0 ] ); 
    % semilogy( sol.it.mtot, sol.metrics.meas_OUT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 

    hold off
    grid on
    title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    
    %legend
    % legend({'total', 'IN random subset',  'OUT random subset'})

%     subplot( 2, 1, 2 ); 
%     hold on
% 
%     semilogy( sol.it.mtot, sol.metrics.exwv_change( : ) , '-o', 'Linewidth', 2 )
% 
% 
%     % semilogy( sol.it.mtot, sol.metrics.exwv_SP, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
%     % semilogy( sol.it.mtot, sol.metrics.exwv_SP_IN, '-o', 'Linewidth', 2, 'Color', [0.0, 0.8, 0.0 ] ); 
%     % semilogy( sol.it.mtot, sol.metrics.exwv_SP_OUT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
% 
%     hold off
%     % title('$ \sum_s || \phi_s - P( \mathbf{r} )  T( \mathbf{r} - \mathbf{r}_s ) ||_F $','Interpreter','latex');
%     grid on
%     title('$ \frac{1}{N_p} \frac{1}{N_s} \sum_s \sum_p \left \Vert \psi_{sp} -  \phi_p \odot T_s \right\Vert^2_F $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     % legend({'total', 'IN random subset',  'OUT random subset'})


    % export_fig( num2str( sol.it.exwv, 'meas_metric-%d.jpg' ), '-r90.0' )
    export_fig( 'meas_metric.jpg', '-r90.0' )

    close all;

    %======================
    % PROBE SCALING METRICS
    %======================
    
%     [ scpm ] = compute_scpm_photonocc( sol.probe.phi );
% 
%     sol.metrics.scpm_fro2TOT(      sol.it.metr ) = scpm.fro2TOT;
%     sol.metrics.scpm_fro2dominant( sol.it.metr ) = scpm.fro2( end );
%     sol.metrics.scpm_fro2others(   sol.it.metr ) = scpm.fro2TOT - scpm.fro2( end );
%     sol.metrics.scpmocc_dominant(  sol.it.metr ) = sol.metrics.scpm_fro2dominant( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
%     sol.metrics.scpmocc_others(    sol.it.metr ) = sol.metrics.scpm_fro2others( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
% 
%     figure( 666 ); 
%     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )  
%     subplot( 4, 1, 1 ); 
%     semilogy( sol.it.mtot, sol.metrics.scpm_fro2TOT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.6, 0.0 ] ); 
%     grid on
%     title('Total Fro norm of probe modes');
%     subplot( 4, 1, 2 ); 
%     semilogy( sol.it.mtot, sol.metrics.scpm_fro2dominant, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
%     grid on
%     title('Fro norm of dominant scpm');
%     subplot( 4, 1, 3 ); 
%     semilogy( sol.it.mtot, sol.metrics.scpm_fro2others, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
%     grid on
%     title('Total Fro norm of other scpm');
%     subplot( 4, 1, 4 ); 
%     hold on
%     semilogy( sol.it.mtot, sol.metrics.scpmocc_dominant, '-o', 'Linewidth', 2, 'Color', [0.5, 0.0, 0.8 ] ); 
%     semilogy( sol.it.mtot, sol.metrics.scpmocc_others, '-o', 'Linewidth', 2, 'Color', [0.2, 0.6, 0.4 ] ); 
%     hold off
%     title('Occupancy of dominant vs others');
%     grid on
%     export_fig( 'metrics_probe_scaling.jpg', '-r90.0' )
% 
%     close all;

    %=====================================================================
    % update the counter that keeps track of metric computation occurances
    %=====================================================================
    
    sol.it.metr = sol.it.metr + 1;   

end
