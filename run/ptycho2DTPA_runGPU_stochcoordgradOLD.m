function [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs )

    %=======================================================
    % form initial exitwaves from current sample/probe modes
    %=======================================================

    sol.sample_sposview_indices = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );
    
%     TF = sol.sample.TF( : );
%     TF = reshape( TF( sol.sample_sposview_indices ), [ sol.sz.sz, 1, size( sol.sample_sposview_indices, 2 ) ]);
    
    TF = reshape( sol.sample.TF( sol.sample_sposview_indices ), [ sol.sz.sz, 1, size( sol.sample_sposview_indices, 2 ) ]);
    
    sol.phi = TF .* sol.probe.P;

    clear( 'TF' )

    %=========================================================
    % perform data movement from CPU main memory to GPU memory
    %=========================================================

    if sol.use_gpu == true

        [ sol.GPU ] = CPUmem2GPUmem( sol, expt );

    else

        [ sol.GPU ] = CPUmem( sol, expt );

    end
        
    %===============================================
    % perform exitwave, sample, probe updates on GPU
    %===============================================

%         exwv_diff = zeros( 1, N_exwvup, 'single' );
    
    for kk = 1 : N_epochs
        
        tic
        
        %===============================
        % Feedback for iteration counter
        %===============================

        S2 = num2str( [ sol.it.epoch, 0, 0, 0, 0 ], 'Epoch = %d, Iteration = %d / %d, Update Repeat = %d / %d' );

        fprintf( [ S2, '\n'])

        %================
        % Exitwave Update 
        %================

%             GPU_phiOLD = gather( sol.GPU.phi );   % USE FOR METRICS, DON'T RECOMPUTE LATER !!

        tic    
        
        %========

%         sol.GPU.phi = ER_GPU_arrays_hadamard( sol.GPU.probe,     ...                      % ERvec
%                                               sol.GPU.TF( : ),   ...
%                                               sol.GPU.ind,       ...
%                                               sol.GPU.sz,        ...
%                                               sol.GPU.Nspos,     ...
%                                               sol.GPU.sqrt_rc,   ...
%                                               sol.GPU.meas_D,    ...
%                                               sol.GPU.meas_Deq0, ...
%                                               sol.GPU.measLPF );

        %========
        
        sol.GPU.phi = RAAR_GPU_arrays_hadamard( sol.GPU.phi,         ...                % RAARvec
                                                sol.GPU.probe,       ...  
                                                sol.GPU.TF( : ),     ...
                                                sol.GPU.ind,         ...
                                                sol.GPU.sz,          ...
                                                sol.GPU.Nspos,       ...
                                                sol.GPU.sqrt_rc,     ...
                                                sol.GPU.meas_D,      ...
                                                sol.GPU.meas_Deq0,   ...
                                                sol.GPU.measLPF,     ...
                                                sol.GPU.RAAR_beta );
                                                   
        %========
        
%         sol.GPU.phi = RAAR_GPU_arrays_hadamard_v2( sol.GPU.phi,         ...              % avgDMERvec
%                                                    sol.GPU.probe,       ...  
%                                                    sol.GPU.TF( : ),     ...
%                                                    sol.GPU.ind,         ...
%                                                    sol.GPU.sz,          ...
%                                                    sol.GPU.Nspos,       ...
%                                                    sol.GPU.sqrt_rc,     ...
%                                                    sol.GPU.meas_D,      ...
%                                                    sol.GPU.meas_Deq0,   ...
%                                                    sol.GPU.measLPF,     ...
%                                                    sol.GPU.RAAR_beta );

        %=======
        
        sol.dt_exwv_up( sol.it.epoch ) = toc;
        
%             exwv_diff( gg ) = norm( gather( sol.GPU.phi( : )) - GPU_phiOLD( : ));

        %============================================================================
        % Now, update the sample transmission and SCPMs using these exit wave updates
        %============================================================================

        %%%%%%%%%%%%%%%
        % Sample Update
        %%%%%%%%%%%%%%%
        
        tic
        
        % For stochastic gradient descent, scramble the sequential update order
        update_order = uint32( gpuArray.randperm( sol.GPU.Nspos ));
        
        if mod( sol.it.epoch, sol.it.sample_update ) == 0

            sol.GPU.TF = ePIEupdate_sample( sol.GPU.phi,   ...       
                                            sol.GPU.probe, ...
                                            sol.GPU.TF,    ...
                                            sol.GPU.vs_r,  ...
                                            sol.GPU.vs_c,  ...
                                            sol.GPU.rs,    ...
                                            update_order,  ...  
                                            'px' );       

            %====================
            % Sparse Sample Edges
            %====================

            if mod( sol.it.epoch, sol.it.sample_sparsity ) == 0

                %==========================
                % Shrinkwrap sample support
                %==========================

    %             TF  = reshape( sol.GPU.TFv, sol.GPU.samsz );
    %             TF  = lpf_gauss( TF, 0.70 * sol.GPU.samsz );
    %             TF  = sparseFDxy_update_2Dsample( TF, sol.GPU.sparseGPU );      
    %             sol.GPU.TFv = TF( : );   
    %             clear( 'TF' )

                %===========================
                % Sparsity using TV variants
                %===========================

                [ sol.GPU.TF ] = sparseFDxy_update_2Dsample( sol.GPU.TF, sol.GPU.sparseGPU );

            end

            %======================================================
            % Sample inequality constraints (projection operations)
            %======================================================

            if mod( sol.it.epoch, sol.it.sample_mag_phs_ineq ) == 0

                sol.GPU.TF = modulus_limits_project( sol.GPU.TF, sol.GPU.abs_TF_lim );

    %             sol.GPU.TF = modulus_limits_project( sol.GPU.TF, [ 0, 1 ] );
    %             sol.GPU.TF = modulus_limits_scale( sol.GPU.TF, sol.GPU.abs_TF_lim );

%                     sol.GPU.TF = phase_limits_project( sol.GPU.TF, sol.GPU.phs_TF_lim );
%                     sol.GPU.TF = phase_limits_scale( sol.GPU.TF, sol.GPU.phs_TF_lim );

            end

            %======================================================
            % Sample mag and phase correlation ( weighted average )
            %======================================================

%             if 0 %mod( sol.it.epoch, 50 ) == 0
% 
%                 abs_TF = abs( sol.GPU.TFv );
%                 abs_TF = abs_TF - min( abs_TF( : ));
%                 abs_TF = abs_TF / max( abs_TF( : ));
% 
%                 phs_TF = angle( sol.GPU.TFv );
%                 phs_TF = phs_TF - min( phs_TF( : ));
%                 phs_TF = phs_TF / max( phs_TF( : ));
% 
%     %             as = 0.75;
%     %             ps = 0.9;
%     %             sol.GPU.TF = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
% 
%                 as = 0.5;
%                 sol.GPU.TF = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.TF ));
% 
%             end  

        end
        
        sol.dt_sample_up( sol.it.epoch ) = toc;
        
        sol.dt_sample_up( sol.it.epoch )
        
        %%%%%%%%%%%%%%
        % Probe Update
        %%%%%%%%%%%%%%
        
        tic
        
        % For stochastic gradient descent, scramble the sequential update order
        update_order = uint32( gpuArray.randperm( sol.GPU.Nspos ));
        
        if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch > sol.it.probe_start )

            sol.GPU.probe = ePIEupdate_probemodes( sol.GPU.phi,      ...
                                                   sol.GPU.probe,    ...
                                                   sol.GPU.TF,       ...
                                                   sol.GPU.vs_r,     ...
                                                   sol.GPU.vs_c,     ...
                                                   sol.GPU.rs,       ...
                                                   update_order,     ...  
                                                   'px' );
            
            %==============
            % Probe Support
            %==============

            if ( mod( sol.it.epoch, sol.it.probe_support ) == 0 ) 

                %=========================
                % Fixed predefined support
                %=========================

%                 sol.GPU.probe = sol.GPU.probe .* sol.GPU.probe_support; 

                %=============================================
                % Shrinkwrap support from probe mode intensity
                %=============================================

                tmp0          = sum( abs( sol.GPU.probe ) .^ 2, 3 );
                [ ~, supp ]   = shrinkwrap( tmp0, sol.GPU.swparams ); 
                sol.GPU.probe = sol.GPU.probe .* supp;

                %==================================
                % Support using dominant probe mode
                %==================================

%                     [ sol.GPU.probe( :, :, 3 ), supp ] = shrinkwrap( sol.GPU.probe( :, :, 3 ), swparams ); 
%                     sol.GPU.probe = sol.GPU.probe .* supp;

            end

            %===================
            % Max abs constraint
            %===================

%             if 0 %( mod( sol.it.epoch, 1 ) == 0 ) && ~isempty( sol.GPU.scpmmax )
% 
% %                     tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
%                 tmp2 = sol.GPU.scpmmax;
%                 tmp0 = ( abs( sol.GPU.probe ) > tmp2 );
%                 tmp1 = not( tmp0 );
% 
%                 sol.GPU.probe = sol.GPU.probe .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.probe )) .* tmp0;
% 
%                 clear( 'tmp0', 'tmp1' )
% 
%             end

            %===========================================================
            % Probe scaling ( # Photons ) and SCPM occupancy constraints
            %===========================================================

            if ( mod( sol.it.epoch, sol.it.probe_scaling ) == 0 ) 

                [ sol.GPU.probe, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.probe, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

            end

        end
        
        sol.dt_probe_up( sol.it.epoch ) = toc;
        
        %========================================================
        % Orthogonalize the SCPMs (spatial coherence probe modes)
        %========================================================

        if ( mod( sol.it.epoch, sol.it.probe_orthog ) == 0 )

            [ sol.GPU.probe ] = orthog_modes_eigendecomp( sol.GPU.probe ); 

        end

        %===========================
        % Fixed Wall time accounting
        %===========================

%         t = toc;
%         sol.epoch_t( sol.it.epoch ) = t;
%         sol.total_t = sol.total_t + t;
% 
%         if sol.total_t > sol.total_t_max
% 
%             [ sol ] = ptycho2DTPA_collectmetrics( sol, expt, sol.spos.batch_N_iteration, kk );
% 
%             [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );
% 
%             sol.it.epoch = sol.it.epoch + 1;  
% 
%             sol.sample.TF = gather( reshape( sol.GPU.TFv, [ sol.sample.sz.sz ]));
%             sol.probe.P   = gather( sol.GPU.probe );
% 
%             sol = rmfield( sol, 'phi' );
%             sol = rmfield( sol, 'GPU' );
% 
%             return
%         end
        
        sol.dt_epoch( sol.it.epoch ) = toc;
        
        %==============================
        % Collect cost function metrics
        %==============================

        if ( mod( sol.it.epoch, sol.it.collect_metrics ) == 0 ) || ( sol.it.epoch == 1e99 )

            [ sol ] = ptycho2DTPA_collectmetrics( sol, expt );

        end

        %=============
        % Plot results
        %=============

        if ( mod( sol.it.epoch, sol.it.print_img_results ) == 0 ) || ( sol.it.epoch == 1e99 )

            %=======================
            % GPU to CPU main memory
            %=======================

            sol.sample.TF = gather( sol.GPU.TF );
            sol.probe.P   = gather( sol.GPU.probe );

            [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );

        end

        %========================
        % Epoch iteration counter
        %========================

        sol.it.epoch = sol.it.epoch + 1;  

    end

    sol = rmfield( sol, 'sample_sposview_indices' );
    sol = rmfield( sol, 'phi' );
    sol = rmfield( sol, 'GPU' );

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    
    %==============
    
    GPU.RAAR_beta = gpuArray( sol.RAAR.beta );
    
    GPU.phi = gpuArray( sol.phi );                                                       

    %==============
    
    GPU.rs    = gpuArray( sol.spos.rs );
    GPU.Nspos = gpuArray( sol.spos.N );
    
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    GPU.sz      = gpuArray( sol.sz.sz );
    GPU.rc      = gpuArray( sol.sz.rc );
    GPU.sqrt_rc = gpuArray( sol.sz.sqrt_rc );

    GPU.samsz = gpuArray( sol.sample.sz.sz );
    GPU.samrc = gpuArray( sol.sample.sz.rc );

    %==============

    GPU.ind        = uint32( sol.sample_sposview_indices ); 
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
    GPU.ind        = gpuArray( GPU.ind ); 
    GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %==============
    
    GPU.meas_D = reshape( expt.meas.D, [ GPU.sz, 1, GPU.Nspos ] );

    GPU.meas_D    = gpuArray( GPU.meas_D );
    GPU.meas_Deq0 = gpuArray( GPU.meas_D == 0 );
    GPU.measLPF   = gpuArray( sol.measLPF );
    
    %==============

    GPU.probe         = gpuArray( sol.probe.P );
    GPU.fro2TOT       = gpuArray( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = gpuArray( sol.probe.scpm.occ );
    GPU.probe_support = gpuArray( sol.probe.support );
    GPU.scpmmax       = gpuArray( sol.probe.scpm.max );
    
    %==============
    
    GPU.swparams.blurx     = gpuArray( sol.swparams.blurx ); 
    GPU.swparams.blury     = gpuArray( sol.swparams.blury );
    GPU.swparams.sparselvl = gpuArray( sol.swparams.sparselvl ); 
    
    %==============
    
    GPU.TF         = gpuArray( sol.sample.TF );
    GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);
    
    GPU.vs_r = gpuArray( sol.sample.vs.r );
    GPU.vs_c = gpuArray( sol.sample.vs.c );
    
    %==============

    GPU.sparseGPU.threshname          = sol.sparse.threshname;
    GPU.sparseGPU.threshtype          = sol.sparse.threshtype;
    GPU.sparseGPU.lvl                 = gpuArray( sol.sparse.lvl );
    GPU.sparseGPU.support             = gpuArray( sol.sparse.support ); 
    GPU.sparseGPU.s2DFDxy.fft_scaling = gpuArray( sol.sparse.s2DFDxy.fft_scaling );
    GPU.sparseGPU.s2DFDxy.qLPF        = gpuArray( sol.sparse.s2DFDxy.qLPF );
    GPU.sparseGPU.s2DFDxy.Dx_fft      = gpuArray( sol.sparse.s2DFDxy.Dx_fft );
    GPU.sparseGPU.s2DFDxy.Dy_fft      = gpuArray( sol.sparse.s2DFDxy.Dy_fft );
    GPU.sparseGPU.s2DFDxy.Dx_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dx_fft_conj );
    GPU.sparseGPU.s2DFDxy.Dy_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dy_fft_conj );
    GPU.sparseGPU.s2DFDxy.sqrt_rc     = gpuArray( sol.sparse.s2DFDxy.sqrt_rc );

    %==============

    % parallel.gpu.GPUArray.zeros(size(data), 'single');

    GPU.spos_conjP_exwv = gpuArray( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
    GPU.spos_abs_P_abs2 = gpuArray( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));

    %====================================================
    % Ensure all relevant information is single precision
    %====================================================
    
    % ???

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )


    %==========================================
    % Move relevant information onto GPU memory
    
    %==============
    
    GPU.RAAR_beta = single( sol.RAAR.beta );
    
    GPU.phi = single( sol.phi );                                                       

    %==============
    
    GPU.rs    = single( sol.spos.rs );
    GPU.Nspos = single( sol.spos.N );
    
    GPU.Nscpm = single( sol.probe.scpm.N );

    GPU.sz      = single( sol.sz.sz );
    GPU.rc      = single( sol.sz.rc );
    GPU.sqrt_rc = single( sol.sz.sqrt_rc );

    GPU.samsz = single( sol.sample.sz.sz );
    GPU.samrc = single( sol.sample.sz.rc );

    %==============

    GPU.ind        = uint32( sol.sample_sposview_indices ); 
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
%     GPU.ind        = gpuArray( GPU.ind ); 
%     GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %==============
    
    GPU.meas_D = reshape( expt.meas.D, [ GPU.sz, 1, GPU.Nspos ] );

    GPU.meas_D    = single( GPU.meas_D );
    GPU.meas_Deq0 = single( GPU.meas_D == 0 );
    GPU.measLPF   = single( sol.measLPF );
    
    %==============

    GPU.probe         = single( sol.probe.P );
    GPU.fro2TOT       = single( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = single( sol.probe.scpm.occ );
    GPU.probe_support = single( sol.probe.support );
    GPU.scpmmax       = single( sol.probe.scpm.max );
    
    %==============
    
    GPU.swparams.blurx     = single( sol.swparams.blurx ); 
    GPU.swparams.blury     = single( sol.swparams.blury );
    GPU.swparams.sparselvl = single( sol.swparams.sparselvl ); 
    
    %==============
    
    GPU.TF         = single( sol.sample.TF );
    GPU.abs_TF_lim = single( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = single( [ sol.sample.phsL, sol.sample.phsH ]);
    
    GPU.vs_r = single( sol.sample.vs.r );
    GPU.vs_c = single( sol.sample.vs.c );
    
    %==============

    GPU.sparseGPU.threshname          = sol.sparse.threshname;
    GPU.sparseGPU.threshtype          = sol.sparse.threshtype;
    GPU.sparseGPU.lvl                 = single( sol.sparse.lvl );
    GPU.sparseGPU.support             = single( sol.sparse.support ); 
    GPU.sparseGPU.s2DFDxy.fft_scaling = single( sol.sparse.s2DFDxy.fft_scaling );
    GPU.sparseGPU.s2DFDxy.qLPF        = single( sol.sparse.s2DFDxy.qLPF );
    GPU.sparseGPU.s2DFDxy.Dx_fft      = single( sol.sparse.s2DFDxy.Dx_fft );
    GPU.sparseGPU.s2DFDxy.Dy_fft      = single( sol.sparse.s2DFDxy.Dy_fft );
    GPU.sparseGPU.s2DFDxy.Dx_fft_conj = single( sol.sparse.s2DFDxy.Dx_fft_conj );
    GPU.sparseGPU.s2DFDxy.Dy_fft_conj = single( sol.sparse.s2DFDxy.Dy_fft_conj );
    GPU.sparseGPU.s2DFDxy.sqrt_rc     = single( sol.sparse.s2DFDxy.sqrt_rc );

    %==============

    % parallel.gpu.single.zeros(size(data), 'single');

    GPU.spos_conjP_exwv = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
    GPU.spos_abs_P_abs2 = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));

    %====================================================
    % Ensure all relevant information is single precision
    %====================================================
    
    % ???
        
end

%====================================================================================================================================================

function [ sol ] = ptycho2DTPA_collectmetrics( sol, expt )

    %=================================
    % get ready for array broadcasting
    %=================================

    meas_D    = reshape( expt.meas.D, [ expt.sz.sz, 1, expt.spos.N ] );
    meas_Deq0 = not( meas_D == 0 );

    %=========================================
    % compute exitwaves for all scan positions
    %=========================================

    sample_sposview_indices = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );

    %========

    TFv  = sol.sample.TF( : );
    TF   = reshape( TFv( sample_sposview_indices ), [ sol.sz.sz, 1, sol.spos.N ]);
    tmp0 = TF .* sol.probe.P;
    tmp1 = fft( fft( fftshift( fftshift( tmp0, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;

    %===========================================================================
    % standard Gaussian noise metric with constant (ignored) stdev at all pixels
    %===========================================================================

    meas_residual = meas_Deq0 .* sqrt( sum( abs( tmp1 ) .^ 2, 3 )) - meas_D;
%     meas_residual = squeeze( sqrt( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 )));
    meas_residual = squeeze( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 ));
    
    clear( 'tmp1', 'TF' )

    %============================
    % RAAR exitwave change metric
    %============================

    phi = RAAR_GPU_arrays_hadamard_v2( tmp0,                    ...
                                       sol.probe.P,             ...
                                       TFv,                     ...
                                       sample_sposview_indices, ...
                                       sol.sz.sz,               ...
                                       sol.spos.N,              ...
                                       sol.sz.sqrt_rc,          ...
                                       meas_D,                  ...
                                       meas_Deq0,               ...
                                       sol.measLPF,             ...
                                       sol.RAAR.beta );


    raar_exwv_change = abs( phi - tmp0 ) .^ 2;
    raar_exwv_change = raar_exwv_change( : );
    
    clear( 'tmp0', 'TFv', 'meas_D', 'meas_Deq0', 'spos' )

    %================
    % collect metrics
    %================

    sol.metrics.exwv_change( sol.it.metr ) = sum( raar_exwv_change ) / length( raar_exwv_change );
    sol.metrics.meas_all( sol.it.metr )    = sum( meas_residual )    / length( meas_residual );
    
    %========

    % fprintf( [ '\n\n', num2str( [ kk, Nit, sol.it.epoch, sol.metrics.meas( sol.it.metr ), sol.metrics.exwv_change( sol.it.metr )], ...
    %             'iteration = %d / %d, iter total = %d, meas metric = %.2f, exit wave sample probe difference = %.2f' ), '\n\n' ]);

    sol.it.mtot( sol.it.metr ) = sol.it.epoch;

    %========

    figure( 666 ); 
    set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     

    subplot( 2, 1, 1 ); 

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

    subplot( 2, 1, 2 ); 
    hold on

    semilogy( sol.it.mtot, sol.metrics.exwv_change( : ) , '-o', 'Linewidth', 2 )


    % semilogy( sol.it.mtot, sol.metrics.exwv_SP, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
    % semilogy( sol.it.mtot, sol.metrics.exwv_SP_IN, '-o', 'Linewidth', 2, 'Color', [0.0, 0.8, 0.0 ] ); 
    % semilogy( sol.it.mtot, sol.metrics.exwv_SP_OUT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 

    hold off
    % title('$ \sum_s || \phi_s - P( \mathbf{r} )  T( \mathbf{r} - \mathbf{r}_s ) ||_F $','Interpreter','latex');
    grid on
    title('$ \frac{1}{N_p} \frac{1}{N_s} \sum_s \sum_p \left \Vert \psi_{sp} -  \phi_p \odot T_s \right\Vert^2_F $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    % legend({'total', 'IN random subset',  'OUT random subset'})


    % export_fig( num2str( sol.it.exwv, 'meas_metric-%d.jpg' ), '-r90.0' )
    export_fig( 'metrics_LSmeasPT_LSphiPT.jpg', '-r90.0' )

    close all;

    %========

    [ scpm ] = compute_scpm_photonocc( sol.probe.P );

    sol.metrics.scpm_fro2TOT(      sol.it.metr ) = scpm.fro2TOT;
    sol.metrics.scpm_fro2dominant( sol.it.metr ) = scpm.fro2( end );
    sol.metrics.scpm_fro2others(   sol.it.metr ) = scpm.fro2TOT - scpm.fro2( end );
    sol.metrics.scpmocc_dominant(  sol.it.metr ) = sol.metrics.scpm_fro2dominant( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
    sol.metrics.scpmocc_others(    sol.it.metr ) = sol.metrics.scpm_fro2others( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );

    figure( 666 ); 
    set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )  
    subplot( 4, 1, 1 ); 
    semilogy( sol.it.mtot, sol.metrics.scpm_fro2TOT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.6, 0.0 ] ); 
    grid on
    title('Total Fro norm of probe modes');
    subplot( 4, 1, 2 ); 
    semilogy( sol.it.mtot, sol.metrics.scpm_fro2dominant, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
    grid on
    title('Fro norm of dominant scpm');
    subplot( 4, 1, 3 ); 
    semilogy( sol.it.mtot, sol.metrics.scpm_fro2others, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
    grid on
    title('Total Fro norm of other scpm');
    subplot( 4, 1, 4 ); 
    hold on
    semilogy( sol.it.mtot, sol.metrics.scpmocc_dominant, '-o', 'Linewidth', 2, 'Color', [0.5, 0.0, 0.8 ] ); 
    semilogy( sol.it.mtot, sol.metrics.scpmocc_others, '-o', 'Linewidth', 2, 'Color', [0.2, 0.6, 0.4 ] ); 
    hold off
    title('Occupancy of dominant vs others');
    grid on
    export_fig( 'metrics_probe_scaling.jpg', '-r90.0' )

    close all;

    %=====================================================================
    % update the counter that keeps track of metric computation occurances
    %=====================================================================
    
    sol.it.metr = sol.it.metr + 1;   

end

%====================================================================================================================================================

function [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt )

        close all;
    
        pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sample.sz.c);
        pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
        pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
        
        pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sample.sz.r);
        pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
        pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
        
        h1 = figure();  
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
%         ax1 = subaxis(1,2,1,'MR',0.1, 'ML',0.1); 
        ax1 = subplot(131);
        imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.sample.TF )); 
%         imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.sample.TF ), [ 0, 1 ]); 
%         imagesc( pltopts.xaxis, pltopts.yaxis,  log10(1 + abs(sol.sample.TF))); 
        %daspect([1 1 1]); 
        axis square
        colorbar
        colormap( ax1, expt.cm.blj )
%         colormap gray; 
        grid on; 
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('abs sample')
        
%         ax2 = subaxis(1,2,2,'MR',0.1, 'ML',0.1); 
        ax2 = subplot(132);
        imagesc( pltopts.xaxis, pltopts.yaxis, angle( sol.sample.TF ), [ -pi, pi ] ); 
%         imagesc( pltopts.xaxis, pltopts.yaxis, angle( sol.sample.TF ), [ sol.sample.phsL, sol.sample.phsH ] ); 
        %daspect([1 1 1]); 
        axis square
        colorbar
%         colormap( ax2, expt.cm.blj )
        colormap( ax2, expt.cm.hsvD )
%         colormap hsv; 
        grid on; 
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('phase sample')
        
        
%         subaxis(1,2,2,'MR',0.1, 'ML',0.1); 
        subplot(133)
        imagescHSV( sol.sample.TF, pltopts  ); 
%         imagescHSV( log10(1 + abs( sol.sample.TF )) .* exp( 1i * angle( sol.sample.TF )), pltopts );
        %daspect([1 1 1]); 
        axis square
        grid on;
        title('HSV ( V = mag, H = phs ) sample')
        
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 1.0 )
%         export_fig( num2str( sol.it.exwv, 'sample_%d.jpg' ), '-r120.0' )
        export_fig( num2str( sol.it.epoch, 'sample_%d.jpg' ), '-r120.0' )
        close all;
     
        %==========================================================================================
        
        [ scpm ] = compute_scpm_photonocc( sol.probe.P );

        close all;
        pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sz.c);
        pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
        pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
        
        pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sz.r);
        pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
        pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
        
        h1 = figure();        
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        for pp = 1 : sol.probe.scpm.N

%             subaxis(2,sol.probe.scpm.N,pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
            subplot(2, sol.probe.scpm.N, pp )   
            imagescHSV(sol.probe.P(:,:,pp), pltopts); 
%             imagescHSV( log10( 1 + 15^-1 * abs( sol.probe.P( :, :, pp ))) .* exp( 1i * angle( sol.probe.P( :, :, pp ))), pltopts); 
            %axis square
            daspect([1 1 1]); 
            %axis off
%             title(num2str( [ sol.probe.scpm.occ( pp ), sol.probe.scpm.fro2TOT ], 'HSV ( V = mag, H = phs ), occupancy = %.4f, fro2TOT = %.4f'))
            title( { 'HSV ( V = mag, H = phs )', num2str(  scpm.occ( pp ), 'occupancy = %.4f' ), ...
                                                 num2str(  scpm.fro2TOT, 'fro2TOT = %.4f' ) })
            grid on;
            set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
      
        end
        
        for pp = 1 : sol.probe.scpm.N

%             subaxis(2,sol.probe.scpm.N,sol.probe.scpm.N +pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
            ax( pp ) = subplot(2, sol.probe.scpm.N, sol.probe.scpm.N + pp );
            imagesc( pltopts.xaxis, pltopts.yaxis, abs(sol.probe.P(:,:,pp))); 
%             imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 15^-1 * abs( sol.probe.P( :, :, pp ))) ); 
            %axis square
            daspect([1 1 1]);
%             colormap gray; 
            colormap( ax( pp ), expt.cm.blj )
            colorbar
            grid on;
            set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
            title('abs probe')
        end

%         export_fig( num2str( sol.it.exwv, 'probe_%d.jpg' ), '-r90.0' )
        export_fig( num2str( sol.it.epoch, 'probe_%d.jpg' ), '-r90.0' )
        close all;

        
        

        absPmodes2 = sqrt( sum( abs( sol.probe.P ) .^ 2, 3 ));
        
        [ phi, ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs.r, sol.sample.vs.c, sol.spos.rs( round( 0.5 * expt.spos.N ), : ), sol.spos.shifttype );
%         phi = fftshift( fft2( fftshift( phi ))) / sqrt( numel( phi ));
        V = fft( fftshift( fft( fftshift( phi, 1 ), [], 1 ), 2 ), [], 2 ) / sqrt( numel( phi ));
        V = fftshift( sqrt( sum( abs( V ) .^ 2, 3 )));

        h1 = figure();        
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        a1 = subplot(121);
        imagesc( pltopts.xaxis, pltopts.yaxis, absPmodes2 ); 
%         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
        daspect([1 1 1]);
        colormap( a1, expt.cm.blj ); 
        colorbar
        grid on;
        set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('abs probe')

        a2 = subplot(122);
        imagesc( log10( 1 + abs( V )));
%         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
        daspect([1 1 1]);
        colormap( a2, expt.cm.blj ); 
        colorbar
        grid on;
        set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('abs^2 fft2 of typical exit wave')
        
        
%         export_fig( num2str( sol.it.exwv, 'absPmodes2_%d.jpg' ), '-r90.0' )
        export_fig( num2str( sol.it.epoch, 'absPmodes2_%d.jpg' ), '-r90.0' )
        close all;
        
        
end

%====================================================================================================================================================

