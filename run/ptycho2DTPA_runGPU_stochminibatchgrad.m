function [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs, N_repeat )

    %=================================================
    % Compute initial exitwaves for all scan positions   % NO!!!! CAN'T DO THIS IF 10^6 SPOS, IS MAIN MEMORY LIMITED
    %=================================================

    % !!!! COMPUTE THIS IN BATCHES
    
    % compute_exitwaves_minibatch( )
    
%     sol.sample_sposview_indices = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );
%             
%     sol.psi = sol.probe.phi .* reshape( sol.sample.T( sol.sample_sposview_indices ), [ sol.sz.sz, 1, size( sol.sample_sposview_indices, 2 ) ]);
    
    %================================================================
    % Send to GPU parameters that will remain constant for all epochs
    %================================================================
    
    sol.GPU.rPIE_alpha = gpuArray( sol.rPIE_alpha );
    sol.GPU.RAAR_beta  = gpuArray( sol.RAAR_beta );

    sol.GPU.sparseGPU.threshname          = sol.sparse.threshname;
    sol.GPU.sparseGPU.threshtype          = sol.sparse.threshtype;
    sol.GPU.sparseGPU.lvl                 = gpuArray( sol.sparse.lvl );
    sol.GPU.sparseGPU.support             = gpuArray( sol.sparse.support ); 
    sol.GPU.sparseGPU.s2DFDxy.fft_scaling = gpuArray( sol.sparse.s2DFDxy.fft_scaling );
    sol.GPU.sparseGPU.s2DFDxy.qLPF        = gpuArray( sol.sparse.s2DFDxy.qLPF );
    sol.GPU.sparseGPU.s2DFDxy.Dx_fft      = gpuArray( sol.sparse.s2DFDxy.Dx_fft );
    sol.GPU.sparseGPU.s2DFDxy.Dy_fft      = gpuArray( sol.sparse.s2DFDxy.Dy_fft );
    sol.GPU.sparseGPU.s2DFDxy.Dx_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dx_fft_conj );
    sol.GPU.sparseGPU.s2DFDxy.Dy_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dy_fft_conj );
    sol.GPU.sparseGPU.s2DFDxy.sqrt_rc     = gpuArray( sol.sparse.s2DFDxy.sqrt_rc );    
        
    sol.GPU.swparams.blurx     = gpuArray( sol.swparams.blurx ); 
    sol.GPU.swparams.blury     = gpuArray( sol.swparams.blury );
    sol.GPU.swparams.sparselvl = gpuArray( sol.swparams.sparselvl ); 
    
    sol.GPU.sz      = gpuArray( sol.sz.sz );
    sol.GPU.rc      = gpuArray( sol.sz.rc );
    sol.GPU.sqrt_rc = gpuArray( sol.sz.sqrt_rc );

    sol.GPU.samsz = gpuArray( sol.sample.sz.sz );
    sol.GPU.samrc = gpuArray( sol.sample.sz.rc );

    sol.GPU.measLPF = reshape( sol.measLPF, [ sol.GPU.sz, 1, 1 ] );
    sol.GPU.measLPF = gpuArray( sol.GPU.measLPF );
    
    sol.GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    sol.GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);
    
    sol.GPU.vs_r = gpuArray( sol.sample.vs.r );
    sol.GPU.vs_c = gpuArray( sol.sample.vs.c );
    
    sol.GPU.T_vec = gpuArray( sol.sample.T( : ));
    
    sol.GPU.phi           = gpuArray( sol.probe.phi );
    sol.GPU.fro2TOT       = gpuArray( sol.probe.scpm.fro2TOT );
    sol.GPU.scpmocc       = gpuArray( sol.probe.scpm.occ );
    sol.GPU.probe_support = gpuArray( sol.probe.support );
    sol.GPU.scpmmax       = gpuArray( sol.probe.scpm.max );

    %==========
    % Main Loop
    %==========

    sol.spos.batch_N = single( round( sol.spos.rand_spos_subset_pct * sol.spos.N ));

    for kk = 1 : N_epochs
        
        %=================================================================================================
        % split the scan position indices up so that we update all positions once to complete a full epoch
        %=================================================================================================
    
        % scamble the index order for all scan position 
        scramble_all_spos_indx = single( randperm( sol.spos.N ));

        % determine number of iterations necessary to complete full epoch
        sol.spos.batch_N_iteration = ceil( sol.spos.N / sol.spos.batch_N );

        % initialize cell to contain mini-batch indices
        sol.spos.batch_nonrepeat_indx = {};
        
        % populate the mini-batches with relevant indices
        if sol.spos.batch_N_iteration == 1
            
            sol.spos.batch_nonrepeat_indx{ 1 } = scramble_all_spos_indx;
            
        else

            for ii = 1 : ( -1 + sol.spos.batch_N_iteration )

                sol.spos.batch_nonrepeat_indx{ ii } = scramble_all_spos_indx( (( ii - 1 ) * sol.spos.batch_N + 1 ) : ( ii * sol.spos.batch_N ));

            end

            ii = sol.spos.batch_N_iteration;

            sol.spos.batch_nonrepeat_indx{ ii } = scramble_all_spos_indx( (( ii - 1 ) * sol.spos.batch_N + 1 ) : end  );

        end

        %=================================
        % now, loop over these minibatches
        %=================================

        for uu = 1 : sol.spos.batch_N_iteration

            %=========================================================
            % perform data movement from CPU main memory to GPU memory
            %=========================================================
            
            sol.spos.batch_indx = sol.spos.batch_nonrepeat_indx{ uu };

            if sol.use_gpu == true

                [ sol.GPU ] = CPUmem2GPUmem( sol, expt, sol.GPU );

            else

                [ sol.GPU ] = CPUmem( sol, expt, sol.GPU );

            end

            %===============================================
            % perform exitwave, sample, probe updates on GPU
            %===============================================

    %         exwv_diff = zeros( 1, N_repeat, 'single' );
    
            for gg = 1 : N_repeat

                %===============================
                % Feedback for iteration counter
                %===============================

                S2 = num2str( [ kk, N_epochs, sol.it.epoch, uu, sol.spos.batch_N_iteration, gg, N_repeat ], 'Local Epoch = %d / %d, Total Epochs = %d, Iteration = %d / %d, Update Repeat = %d / %d' );

                fprintf( [ S2, '\n'])
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                   Exitwave Update
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tic
                
    %             GPU_psiOLD = gather( sol.GPU.psi );
    
                %========
                
                sol.GPU.psi = ER_GPU_arrays_hadamard( sol.GPU.phi,       ...            % ERvec_2DTPA
                                                      sol.GPU.T_vec,     ...
                                                      sol.GPU.ind,       ...
                                                      sol.GPU.sz,        ...
                                                      sol.GPU.Nspos,     ...
                                                      sol.GPU.sqrt_rc,   ...
                                                      sol.GPU.meas_D,    ...
                                                      sol.GPU.meas_Deq0, ...
                                                      sol.GPU.measLPF );
  
                %========

%                 sol.GPU.psi = RAAR_GPU_arrays_hadamard( sol.GPU.psi,         ...                % RAARvec_2DTPA
%                                                         sol.GPU.phi,         ...  
%                                                         sol.GPU.T_vec,       ...
%                                                         sol.GPU.ind,         ...
%                                                         sol.GPU.sz,          ...
%                                                         sol.GPU.Nspos,       ...
%                                                         sol.GPU.sqrt_rc,     ...
%                                                         sol.GPU.meas_D,      ...
%                                                         sol.GPU.meas_Deq0,   ...
%                                                         sol.GPU.measLPF,     ...
%                                                         sol.GPU.RAAR_beta );

                %========

%                 sol.GPU.psi = RAAR_GPU_arrays_hadamard_v2( sol.GPU.psi,         ...          % mRAARvec_2DTPA
%                                                            sol.GPU.phi,         ...  
%                                                            sol.GPU.T_vec,       ...
%                                                            sol.GPU.ind,         ...
%                                                            sol.GPU.sz,          ...
%                                                            sol.GPU.Nspos,       ...
%                                                            sol.GPU.sqrt_rc,     ...
%                                                            sol.GPU.meas_D,      ...
%                                                            sol.GPU.meas_Deq0,   ...
%                                                            sol.GPU.measLPF,     ...
%                                                            sol.GPU.RAAR_beta );

                %=======

    %             exwv_diff( gg ) = norm( gather( sol.GPU.psi( : )) - GPU_psiOLD( : ));
    
                sol.timings.exwv_update( sol.it.epoch ) = toc;
   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                       Sample Update
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                tic
                
                if mod( sol.it.epoch, sol.it.sample_update ) == 0
           
                    [ sol.GPU.T_vec ] = batchgradupdate_2DTPAsample( sol.GPU.psi,        ...
                                                                     sol.GPU.T_vec,      ...
                                                                     sol.GPU.phi,        ...
                                                                     sol.GPU.ind_offset, ...
                                                                     sol.GPU.rc,         ...
                                                                     sol.GPU.Nspos,      ...
                                                                     sol.GPU.rPIE_alpha );            
                    %====================
                    % Sparse Sample Edges
                    %====================

                    if mod( sol.it.epoch, sol.it.sample_sparsity ) == 0

                        %==========================
                        % Shrinkwrap sample support
                        %==========================
                        
            %             TF  = reshape( sol.GPU.T_vec, sol.GPU.samsz );
            %             TF  = lpf_gauss( TF, 0.70 * sol.GPU.samsz );
            %             TF  = sparseFDxy_update_2Dsample( TF, sol.GPU.sparseGPU );      
            %             sol.GPU.T_vec = TF( : );   
            %             clear( 'TF' )

                        %===========================
                        % Sparsity using TV variants
                        %===========================
                        
                        [ sol.GPU.T_vec ] = sparseFDxy_update_2Dsample( reshape( sol.GPU.T_vec, sol.GPU.samsz ), sol.GPU.sparseGPU );

                        sol.GPU.T_vec = sol.GPU.T_vec( : );
 
                    end

                    %======================================================
                    % Sample inequality constraints (projection operations)
                    %======================================================

                    if mod( sol.it.epoch, sol.it.sample_mag_phs_ineq ) == 0

                        sol.GPU.T_vec = modulus_limits_project( sol.GPU.T_vec, sol.GPU.abs_TF_lim );

            %             sol.GPU.T_vec = modulus_limits_project( sol.GPU.T_vec, [ 0, 1 ] );
            %             sol.GPU.T_vec = modulus_limits_scale( sol.GPU.T_vec, sol.GPU.abs_TF_lim );

    %                     sol.GPU.T_vec = phase_limits_project( sol.GPU.T_vec, sol.GPU.phs_TF_lim );
    %                     sol.GPU.T_vec = phase_limits_scale( sol.GPU.T_vec, sol.GPU.phs_TF_lim );

                    end

                    %======================================================
                    % Sample mag and phase correlation ( weighted average )
                    %======================================================

%                     if 0 %mod( sol.it.epoch, 50 ) == 0
% 
%                         abs_TF = abs( sol.GPU.T_vec );
%                         abs_TF = abs_TF - min( abs_TF( : ));
%                         abs_TF = abs_TF / max( abs_TF( : ));
% 
%                         phs_TF = angle( sol.GPU.T_vec );
%                         phs_TF = phs_TF - min( phs_TF( : ));
%                         phs_TF = phs_TF / max( phs_TF( : ));
% 
%             %             as = 0.75;
%             %             ps = 0.9;
%             %             sol.GPU.T_vec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
% 
%                         as = 0.5;
%                         sol.GPU.T_vec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.T_vec ));
% 
%                     end  

                end
                
                sol.timings.sample_update( sol.it.epoch ) = toc;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                       Probe Update
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tic
                
                if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch > sol.it.probe_start )

                    %==================================
                    % DM/ePIE hybrid style probe update  
                    %==================================

                    T_view = reshape( sol.GPU.T_vec( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);

                    abs2_TFview = abs( T_view ) .^ 2;

                    sol.GPU.phi = sol.GPU.phi + ( sum( conj( T_view ) .* sol.GPU.psi - sol.GPU.phi .* abs2_TFview, 4 ) ) ./ ( 1e-7 + sum( abs2_TFview, 4 ) );
                    
                    %==============
                    % Probe Support
                    %==============

                    if ( mod( sol.it.epoch, sol.it.probe_support ) == 0 ) 

                        %=============================================
                        % Shrinkwrap support from probe mode intensity
                        %=============================================
                        
                        tmp0          = sum( abs( sol.GPU.phi ) .^ 2, 3 );
                        [ ~, supp ]   = shrinkwrap( tmp0, sol.GPU.swparams ); 
                        sol.GPU.phi   = sol.GPU.phi .* supp;
                        
                        %=========================
                        % Fixed predefined support
                        %=========================
                        
                        sol.GPU.phi = sol.GPU.phi .* sol.GPU.probe_support; 

                        %==================================
                        % Support using dominant probe mode
                        %==================================
                        
    %                     [ sol.GPU.phi( :, :, 3 ), supp ] = shrinkwrap( sol.GPU.phi( :, :, 3 ), swparams ); 
    %                     sol.GPU.phi = sol.GPU.phi .* supp;

                    end

                    %===================
                    % Max abs constraint
                    %===================

    %                 if 0 %( mod( sol.it.epoch, 1 ) == 0 ) && ~isempty( sol.GPU.scpmmax )
    %                     
    % %                     tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
    %                     tmp2 = sol.GPU.scpmmax;
    %                     tmp0 = ( abs( sol.GPU.phi ) > tmp2 );
    %                     tmp1 = not( tmp0 );
    %                         
    %                     sol.GPU.phi = sol.GPU.phi .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.phi )) .* tmp0;
    %                     
    %                     clear( 'tmp0', 'tmp1' )
    %                     
    %                 end

                    %============================
                    % Probe Scaling ( # Photons )
                    %============================

                    if ( mod( sol.it.epoch, sol.it.probe_scaling ) == 0 ) 

                        [ sol.GPU.phi, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.phi, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

                    end

                end
                
                sol.timings.probe_update( sol.it.epoch ) = toc;
                
            end
        
            %=======================
            % GPU to CPU main memory
            %=======================

%             sol.sample.T  = gather( reshape( sol.GPU.T_vec, [ sol.sample.sz.sz ]));
%             sol.probe.phi = gather( sol.GPU.phi );
            
%             % only need this if using RAAR
%             sol.psi( :, :, :, sol.spos.batch_indx ) = gather( sol.GPU.psi );

        end
                    
%         %============================
%         % Center the probe intensity:
%         %============================
% 
%         if ( mod( sol.it.epoch, sol.it.probe_centering ) == 0 ) 
% 
%             [ com ] = centerofmass( sum( abs( sol.GPU.phi ) .^ 2, 3 ));
%             sol.GPU.phi = circshift( sol.GPU.phi, -1 * round( [ com - 0.5 * sol.GPU.sz - 1, 0] ));
% 
%         end

        %========================================================
        % Orthogonalize the SCPMs (spatial coherence probe modes)
        %========================================================
        
        if ( mod( sol.it.epoch, sol.it.probe_orthog ) == 0 ) || ( sol.it.epoch == 1 )
            
            [ sol.GPU.phi ] = orthog_modes_eigendecomp( sol.GPU.phi ); 
      
        end
   
        %============================================
        % Collect cost function metrics, Plot results
        %============================================
        
        if ( mod( sol.it.epoch, sol.it.metrics_and_plotting ) == 0 ) || ( sol.it.epoch == 1 )

            sol.sample.T  = gather( reshape( sol.GPU.T_vec, [ sol.sample.sz.sz ]));
            sol.probe.phi = gather( sol.GPU.phi );
            
            [ sol ] = ptycho2DTPA_collectmetrics( sol, expt );
            
            ptycho2DTPA_plotresults( sol, expt );

        end
        
        %========================
        % Epoch iteration counter
        %========================

        sol.it.epoch = sol.it.epoch + 1;  
                
    end
    
%     sol = rmfield( sol, 'psi' );
    sol = rmfield( sol, 'GPU' );

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt, GPU )

    GPU.batch_indx = gpuArray( sol.spos.batch_indx );

    GPU.batch_rs = gpuArray( sol.spos.rs( sol.spos.batch_indx, : ) );

    GPU.sample_sposview_indices = get_indices_2Dframes( GPU.batch_rs, sol.GPU.samsz, sol.GPU.vs_r, sol.GPU.vs_c ); 

    GPU.batch_N = gpuArray( length( GPU.batch_indx ) );

    GPU.Nspos = gpuArray( GPU.batch_N );
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    %========

    GPU.ind = uint32( GPU.sample_sposview_indices ); 

    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
    GPU.ind        = gpuArray( GPU.ind ); 
    GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %========
    
    GPU.meas_D = gpuArray( reshape( expt.meas.D( :, :, GPU.batch_indx ), [ GPU.sz, 1, GPU.Nspos ] ));

    GPU.meas_Deq0 = ( GPU.meas_D == 0 );
  
end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )

    GPU = [];

end

%====================================================================================================================================================

% function [ sol ] = ptycho2DTPA_collectmetrics( sol, expt, Nit, kk )
% 
%     %=================================
%     % get ready for array broadcasting
%     %=================================
% 
%     meas_D    = reshape( expt.meas.D, [ expt.sz.sz, 1, expt.spos.N ] );
%     meas_Deq0 = not( meas_D == 0 );
% 
%     %=========================================
%     % compute exitwaves for all scan positions
%     %=========================================
% 
%     spos.batch_indx = 1 : sol.spos.N;
%     spos.batch_rs         = sol.spos.rs( spos.batch_indx, : );
%     spos.frameindx        = get_indices_2Dframes( spos.batch_rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );
% 
%     %========
% 
%     TFv  = sol.sample.T( : );
%     TF   = reshape( TFv( spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N ]);
%     tmp0 = TF .* sol.probe.phi;
%     tmp1 = fft( fft( fftshift( fftshift( tmp0, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;
% 
%     %===========================================================================
%     % standard Gaussian noise metric with constant (ignored) stdev at all pixels
%     %===========================================================================
% 
%     meas_residual = meas_Deq0 .* sqrt( sum( abs( tmp1 ) .^ 2, 3 )) - meas_D;
% %     meas_residual = squeeze( sqrt( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 )));
%     meas_residual = squeeze( sum( sum( abs( meas_residual ) .^ 2, 1 ), 2 ));
%     
%     clear( 'tmp1', 'TF' )
% 
%     %============================
%     % RAAR exitwave change metric
%     %============================
% 
% %     psi = RAAR_GPU_arrays_hadamard( tmp0,            ...
% %                                     sol.probe.phi,     ...
% %                                     TFv,             ...
% %                                     spos.frameindx,  ...
% %                                     sol.sz.sz,       ...
% %                                     sol.spos.N,      ...
% %                                     sol.sz.sqrt_rc,  ...
% %                                     meas_D,          ...
% %                                     meas_Deq0,       ...
% %                                     sol.measLPF,     ...
% %                                     sol.RAAR.beta );
% % 
% % 
% %     raar_exwv_change = abs( psi - tmp0 ) .^ 2;
% % 
% %     clear( 'tmp0', 'TFv', 'meas_D', 'meas_Deq0', 'spos' )
% 
%     %================
%     % collect metrics
%     %================
% 
%     for ii = 1 : length( sol.spos.batch_nonrepeat_indx )
% 
%         meas_residual_batch = meas_residual( sol.spos.batch_nonrepeat_indx{ ii } );
%         meas_residual_batch = meas_residual_batch( : );
% 
%         sol.metrics.meas( sol.it.metr, ii ) = sum( meas_residual_batch ) / length( meas_residual_batch );
% 
% %         raar_exwv_change_batch = raar_exwv_change( sol.spos.batch_nonrepeat_indx{ ii } );
% %         raar_exwv_change_batch = raar_exwv_change_batch( : );
% % 
% %         sol.metrics.exwv_change( sol.it.metr, ii ) = sum( raar_exwv_change_batch ) / length( raar_exwv_change_batch );
% 
%     end
% 
%     sol.metrics.meas_all( sol.it.metr ) = sum( meas_residual ) / length( meas_residual );
%     
%     %========
% 
%     % fprintf( [ '\n\n', num2str( [ kk, Nit, sol.it.epoch, sol.metrics.meas( sol.it.metr ), sol.metrics.exwv_change( sol.it.metr )], ...
%     %             'iteration = %d / %d, iter total = %d, meas metric = %.2f, exit wave sample probe difference = %.2f' ), '\n\n' ]);
% 
%     sol.it.mtot( sol.it.metr ) = sol.it.epoch;
% 
%     %========
% 
%     figure( 666 ); 
%     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     
% 
% %     subplot( 2, 1, 1 ); 
% 
%     hold on
% 
%     for ii = 1 : length( sol.spos.batch_nonrepeat_indx )
% 
%         semilogy( sol.it.mtot, sol.metrics.meas( :, ii ) , '-o', 'Linewidth', 2 )
% 
%     end
%     
%     semilogy( sol.it.mtot, sol.metrics.meas_all, '--', 'Linewidth', 4, 'color', [ 0, 0, 0 ] )
% 
%     % semilogy( sol.it.mtot, sol.metrics.meas, '-o', 'Linewidth', 2, 'Color', [0.8, 0, 0 ] ); 
%     % semilogy( sol.it.mtot, sol.metrics.meas_IN, '-o', 'Linewidth', 2, 'Color', [0.0, 0.8, 0 ] ); 
%     % semilogy( sol.it.mtot, sol.metrics.meas_OUT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
% 
%     hold off
%     grid on
%     title('$ \frac{1}{N_s} \sum_s \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     
%     %legend
%     % legend({'total', 'IN random subset',  'OUT random subset'})
% 
% %     subplot( 2, 1, 2 ); 
% %     hold on
% % 
% % 
% %     for ii = 1 : length( sol.spos.batch_nonrepeat_indx )
% % 
% %         semilogy( sol.it.mtot, sol.metrics.exwv_change( :, ii ) , '-o', 'Linewidth', 2 )
% % 
% %     end
% % 
% %     % semilogy( sol.it.mtot, sol.metrics.exwv_SP, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
% %     % semilogy( sol.it.mtot, sol.metrics.exwv_SP_IN, '-o', 'Linewidth', 2, 'Color', [0.0, 0.8, 0.0 ] ); 
% %     % semilogy( sol.it.mtot, sol.metrics.exwv_SP_OUT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
% % 
% %     hold off
% %     % title('$ \sum_s || \phi_s - P( \mathbf{r} )  T( \mathbf{r} - \mathbf{r}_s ) ||_F $','Interpreter','latex');
% %     grid on
% %     title('$ \frac{1}{N_p} \frac{1}{N_s} \sum_s \sum_p \left \Vert \psi_{sp} -  \phi_p \odot T_s \right\Vert^2_F $', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% %     % legend({'total', 'IN random subset',  'OUT random subset'})
% 
% 
%     % export_fig( num2str( sol.it.exwv, 'meas_metric-%d.jpg' ), '-r90.0' )
%     export_fig( 'metrics_LSmeasPT_LSphiPT.jpg', '-r90.0' )
% 
%     close all;
% 
%     %========
% 
% %     [ scpm ] = compute_scpm_photonocc( sol.probe.phi );
% % 
% %     sol.metrics.scpm_fro2TOT(      sol.it.metr ) = scpm.fro2TOT;
% %     sol.metrics.scpm_fro2dominant( sol.it.metr ) = scpm.fro2( end );
% %     sol.metrics.scpm_fro2others(   sol.it.metr ) = scpm.fro2TOT - scpm.fro2( end );
% %     sol.metrics.scpmocc_dominant(  sol.it.metr ) = sol.metrics.scpm_fro2dominant( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
% %     sol.metrics.scpmocc_others(    sol.it.metr ) = sol.metrics.scpm_fro2others( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
% % 
% %     figure( 666 ); 
% %     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )  
% %     subplot( 4, 1, 1 ); 
% %     semilogy( sol.it.mtot, sol.metrics.scpm_fro2TOT, '-o', 'Linewidth', 2, 'Color', [0.0, 0.6, 0.0 ] ); 
% %     grid on
% %     title('Total Fro norm of probe modes');
% %     subplot( 4, 1, 2 ); 
% %     semilogy( sol.it.mtot, sol.metrics.scpm_fro2dominant, '-o', 'Linewidth', 2, 'Color', [0.8, 0.0, 0.0 ] ); 
% %     grid on
% %     title('Fro norm of dominant scpm');
% %     subplot( 4, 1, 3 ); 
% %     semilogy( sol.it.mtot, sol.metrics.scpm_fro2others, '-o', 'Linewidth', 2, 'Color', [0.0, 0.0, 0.8 ] ); 
% %     grid on
% %     title('Total Fro norm of other scpm');
% %     subplot( 4, 1, 4 ); 
% %     hold on
% %     semilogy( sol.it.mtot, sol.metrics.scpmocc_dominant, '-o', 'Linewidth', 2, 'Color', [0.5, 0.0, 0.8 ] ); 
% %     semilogy( sol.it.mtot, sol.metrics.scpmocc_others, '-o', 'Linewidth', 2, 'Color', [0.2, 0.6, 0.4 ] ); 
% %     hold off
% %     title('Occupancy of dominant vs others');
% %     grid on
% %     export_fig( 'metrics_probe_scaling.jpg', '-r90.0' )
% % 
% %     close all;
% 
%     %=====================================================================
%     % update the counter that keeps track of metric computation occurances
%     %=====================================================================
%     
%     sol.it.metr = sol.it.metr + 1;   
% 
% end

%====================================================================================================================================================

% function [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt )
% 
%   
%         %==============================
%         % probe mode correlation matrix
%         %==============================
% 
%         P = reshape( sol.probe.phi, [ sol.sz.rc, sol.probe.scpm.N ] );
% 
%         corr_matrix_scpm = ctranspose( P ) * P;
%         
%         h1 = figure();  
%         set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%     
%         imagesc( log10( 1 + abs( corr_matrix_scpm )))
%         axis square
%         colorbar
%         colormap( expt.cm.blj )
%         title('Probe Mode Correlation Matrix')
%         
%         export_fig( num2str( sol.it.epoch, 'corr_matrix_scpm_%d.jpg' ), '-r120.0' )
%         close all;
%         
%         %========
% 
%         close all;
%     
%         pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sample.sz.c);
%         pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
%         pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
%         
%         pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sample.sz.r);
%         pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
%         pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
%         
%         h1 = figure();  
%         set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%         
% %         ax1 = subaxis(1,2,1,'MR',0.1, 'ML',0.1); 
%         ax1 = subplot(131);
%         imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.sample.T )); 
% %         imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.sample.T ), [ 0, 1 ]); 
% %         imagesc( pltopts.xaxis, pltopts.yaxis,  log10(1 + abs(sol.sample.T))); 
%         %daspect([1 1 1]); 
%         axis square
%         colorbar
%         colormap( ax1, expt.cm.blj )
% %         colormap gray; 
%         grid on; 
% %         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%         title('abs sample')
%         
% %         ax2 = subaxis(1,2,2,'MR',0.1, 'ML',0.1); 
%         ax2 = subplot(132);
%         imagesc( pltopts.xaxis, pltopts.yaxis, angle( sol.sample.T ), [ -pi, pi ] ); 
% %         imagesc( pltopts.xaxis, pltopts.yaxis, angle( sol.sample.T ), [ sol.sample.phsL, sol.sample.phsH ] ); 
%         %daspect([1 1 1]); 
%         axis square
%         colorbar
% %         colormap( ax2, expt.cm.blj )
%         colormap( ax2, expt.cm.hsvD )
% %         colormap hsv; 
%         grid on; 
% %         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%         title('phase sample')
%         
%         
% %         subaxis(1,2,2,'MR',0.1, 'ML',0.1); 
%         subplot(133)
%         imagescHSV( sol.sample.T, pltopts  ); 
% %         imagescHSV( log10(1 + abs( sol.sample.T )) .* exp( 1i * angle( sol.sample.T )), pltopts );
%         %daspect([1 1 1]); 
%         axis square
%         grid on;
%         title('HSV ( V = mag, H = phs ) sample')
%         
% %         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 1.0 )
% %         export_fig( num2str( sol.it.exwv, 'sample_%d.jpg' ), '-r120.0' )
%         export_fig( num2str( sol.it.epoch, 'sample_%d.jpg' ), '-r120.0' )
%         close all;
%      
%         %==========================================================================================
%         
%         [ scpm ] = compute_scpm_photonocc( sol.probe.phi );
% 
%         close all;
%         pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sz.c);
%         pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
%         pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
%         
%         pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sz.r);
%         pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
%         pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
%         
%         h1 = figure();        
%         set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
%         
%         for pp = 1 : sol.probe.scpm.N
% 
% %             subaxis(2,sol.probe.scpm.N,pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
%             subplot( 2, double( sol.probe.scpm.N ), double( pp ) )   
%             imagescHSV( sol.probe.phi( :, :, pp ), pltopts ); 
% %             imagescHSV( log10( 1 + 15^-1 * abs( sol.probe.phi( :, :, pp ))) .* exp( 1i * angle( sol.probe.phi( :, :, pp ))), pltopts); 
%             %axis square
%             daspect([1 1 1]); 
%             %axis off
% %             title(num2str( [ sol.probe.scpm.occ( pp ), sol.probe.scpm.fro2TOT ], 'HSV ( V = mag, H = phs ), occupancy = %.4f, fro2TOT = %.4f'))
%             title( { 'HSV ( V = mag, H = phs )', num2str(  scpm.occ( pp ), 'occupancy = %.4f' ), ...
%                                                  num2str(  scpm.fro2TOT, 'fro2TOT = %.4f' ) })
%             grid on;
%             set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%       
%         end
%         
%         for pp = 1 : sol.probe.scpm.N
% 
% %             subaxis(2,sol.probe.scpm.N,sol.probe.scpm.N +pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
%             ax( pp ) = subplot( 2, double( sol.probe.scpm.N ), double( sol.probe.scpm.N + pp ) );
%             imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.probe.phi( :, :, pp ))); 
% %             imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 15^-1 * abs( sol.probe.phi( :, :, pp ))) ); 
%             %axis square
%             daspect([1 1 1]);
% %             colormap gray; 
%             colormap( ax( pp ), expt.cm.blj )
%             colorbar
%             grid on;
%             set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%             title('abs probe')
%         end
% 
% %         export_fig( num2str( sol.it.exwv, 'probe_%d.jpg' ), '-r90.0' )
%         export_fig( num2str( sol.it.epoch, 'probe_%d.jpg' ), '-r90.0' )
%         close all;
% 
%         
%         
% 
%         absPmodes2 = sqrt( sum( abs( sol.probe.phi ) .^ 2, 3 ));
%         
%         [ phi, ~ ] = enforce_2DTPAsposview( sol.probe.phi, sol.sample.T, sol.sample.vs.r, sol.sample.vs.c, sol.spos.rs( round( 0.5 * expt.spos.N ), : ), sol.spos.shifttype );
% %         phi = fftshift( fft2( fftshift( phi ))) / sqrt( numel( phi ));
%         V = fft( fftshift( fft( fftshift( phi, 1 ), [], 1 ), 2 ), [], 2 ) / sqrt( numel( phi ));
%         V = fftshift( sqrt( sum( abs( V ) .^ 2, 3 )));
% 
%         h1 = figure();        
%         set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%         
%         a1 = subplot(121);
%         imagesc( pltopts.xaxis, pltopts.yaxis, absPmodes2 ); 
% %         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
%         daspect([1 1 1]);
%         colormap( a1, expt.cm.blj ); 
%         colorbar
%         grid on;
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%         title('abs probe')
% 
%         a2 = subplot(122);
%         imagesc( log10( 1 + abs( V )));
% %         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
%         daspect([1 1 1]);
%         colormap( a2, expt.cm.blj ); 
%         colorbar
%         grid on;
%         set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
%         title('abs^2 fft2 of typical exit wave')
%         
%         
% %         export_fig( num2str( sol.it.exwv, 'absPmodes2_%d.jpg' ), '-r90.0' )
%         export_fig( num2str( sol.it.epoch, 'absPmodes2_%d.jpg' ), '-r90.0' )
%         close all;
%         
%         
% end

%====================================================================================================================================================

