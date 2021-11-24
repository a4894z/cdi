function [ sol, expt ] = ptycho2DTPA_runGPU_stochcoordgrad( sol, expt, N_epochs )

    %=======================================================
    % form initial exitwaves from current sample/probe modes
    %=======================================================

    sol.sample_sposview_indices = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );

%     TF = reshape( sol.sample.T( sol.sample_sposview_indices ), [ sol.sz.sz, 1, size( sol.sample_sposview_indices, 2 ) ]);
%     
%     sol.psi = TF .* sol.probe.phi;
% 
%     clear( 'TF' )

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

    for kk = 1 : N_epochs
        
        start_epoch = tic;

        %============================================================================================================================================
        %                                                      Exitwave Update
        %============================================================================================================================================

        start_exwv = tic;  
        
        %========

        sol.GPU.psi = exitwave_vectorized_update_2DTPA_meas_projection( sol.GPU.phi,     ...                      % ERvec
                                                                        sol.GPU.TF( : ),   ...
                                                                        sol.GPU.ind,       ...
                                                                        sol.GPU.sz,        ...
                                                                        sol.GPU.Nspos,     ...
                                                                        sol.GPU.sqrt_rc,   ...
                                                                        sol.GPU.meas_D,    ...
                                                                        sol.GPU.meas_Deq0, ...
                                                                        sol.GPU.measLPF );
                           
        %========
        
%         sol.GPU.psi = RAAR_GPU_arrays_hadamard_v2( sol.GPU.psi,         ...              % avgDMERvec
%                                                    sol.GPU.phi,       ...  
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
        
        sol.timings.exwv_update( sol.it.epoch ) = toc( start_exwv );

        %============================================================================================================================================
        %                                                         Sample Update
        %============================================================================================================================================
 
        if mod( sol.it.epoch, sol.it.sample_update ) == 0

            start_sample = tic;
      
            %================================================================
            % keep a copy of the current sample for use in probe update below
            %================================================================

            sol.GPU.TF_old = sol.GPU.TF;

            %========
        
            % For stochastic gradient descent, scramble the sequential update order
            update_order = uint32( gpuArray.randperm( sol.GPU.Nspos ));

            sol.GPU.TF = rPIEupdate_blockstoch_2DTPA_sample( sol.GPU.psi,        ...       
                                                             sol.GPU.phi,        ...
                                                             sol.GPU.TF,         ...
                                                             sol.GPU.vs_r,       ...
                                                             sol.GPU.vs_c,       ...
                                                             sol.GPU.rs,         ...
                                                             update_order,       ...  
                                                             sol.GPU.rPIE_alpha, ...
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

            sol.timings.sample_update( sol.it.epoch ) = toc( start_sample );

        end

        %============================================================================================================================================
        %                                                       Exitwave Update
        %============================================================================================================================================
        
%         % Update exitwaves using recently updated sample transfer function
%         start_exwv = tic;  
% 
%         sol.GPU.psi = exitwave_vectorized_update_2DTPA_meas_projection( sol.GPU.phi,     ...                    
%                                                                         sol.GPU.TF( : ),   ...
%                                                                         sol.GPU.ind,       ...
%                                                                         sol.GPU.sz,        ...
%                                                                         sol.GPU.Nspos,     ...
%                                                                         sol.GPU.sqrt_rc,   ...
%                                                                         sol.GPU.meas_D,    ...
%                                                                         sol.GPU.meas_Deq0, ...
%                                                                         sol.GPU.measLPF );
%                                           
%         sol.timings.exwv_update( sol.it.epoch ) = sol.timings.exwv_update( sol.it.epoch ) + toc( start_exwv );

        %============================================================================================================================================
        %                                                                   Probe Update
        %============================================================================================================================================

        if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch > sol.it.probe_start )
            
            start_probe = tic;

            % For stochastic gradient descent, scramble the sequential update order
            update_order = uint32( gpuArray.randperm( sol.GPU.Nspos ));

            %=====================================================================================
            % probe update using new T^{(k+1)} for exitwave update, new T^{(k+1)} for probe update      
            %=====================================================================================
            
%             sol.GPU.phi = rPIEupdate_blockstoch_2DTPA_SCPMprobes( sol.GPU.psi,  ...
%                                                                   sol.GPU.phi,  ...
%                                                                   sol.GPU.TF,   ...
%                                                                   sol.GPU.vs_r, ...
%                                                                   sol.GPU.vs_c, ...
%                                                                   sol.GPU.rs,   ...
%                                                                   update_order, ...  
%                                                                   single( 1.0 ),  ...
%                                                                   'px' );
                                                       
            %=================================================================================
            % probe update using old T^{(k)} for exitwave update, old T^{(k)} for probe update
            %=================================================================================
            
            sol.GPU.phi = rPIEupdate_blockstoch_2DTPA_SCPMprobes( sol.GPU.psi,    ...
                                                                  sol.GPU.phi,    ...
                                                                  sol.GPU.TF_old, ...
                                                                  sol.GPU.vs_r,   ...
                                                                  sol.GPU.vs_c,   ...
                                                                  sol.GPU.rs,     ...
                                                                  update_order,   ...  
                                                                  single( 1.0 ),  ...
                                                                  'px' );                              
                            
            %==============
            % Probe Support
            %==============

            if ( mod( sol.it.epoch, sol.it.probe_support ) == 0 ) 

                %=============================================
                % Shrinkwrap support from probe mode intensity
                %=============================================

                tmp0          = sum( abs( sol.GPU.phi ) .^ 2, 3 );
                [ ~, supp ]   = shrinkwrap( tmp0, sol.GPU.swparams ); 
                sol.GPU.phi = sol.GPU.phi .* supp;

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

%             if 0 %( mod( sol.it.epoch, 1 ) == 0 ) && ~isempty( sol.GPU.scpmmax )
% 
% %                     tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
%                 tmp2 = sol.GPU.scpmmax;
%                 tmp0 = ( abs( sol.GPU.phi ) > tmp2 );
%                 tmp1 = not( tmp0 );
% 
%                 sol.GPU.phi = sol.GPU.phi .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.phi )) .* tmp0;
% 
%                 clear( 'tmp0', 'tmp1' )
% 
%             end

            %===========================================================
            % Probe scaling ( # Photons ) and SCPM occupancy constraints
            %===========================================================

            if ( mod( sol.it.epoch, sol.it.probe_scaling ) == 0 ) 

                [ sol.GPU.phi, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.phi, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

            end
            
%             %============================
%             % Center the probe intensity:
%             %============================
%             
%             if ( mod( sol.it.epoch, sol.it.probe_centering ) == 0 ) 
%                 
%                 [ com ] = centerofmass( sum( abs( sol.GPU.phi ) .^ 2, 3 ));
%                 sol.GPU.phi = circshift( sol.GPU.phi, -1 * round( [ com - 0.5 * sol.GPU.sz - 1, 0] ));
% 
%             end
            
            %========================================================
            % Orthogonalize the SCPMs (spatial coherence probe modes)
            %========================================================

            if ( mod( sol.it.epoch, sol.it.probe_orthog ) == 0 ) || ( sol.it.epoch == 1 )

                [ sol.GPU.phi ] = orthog_modes_eigendecomp( sol.GPU.phi ); 

            end
                    
            sol.timings.probe_update( sol.it.epoch ) = toc( start_probe );

        end

        %=============
        % Epoch timing
        %=============
        
        sol.timings.epoch( sol.it.epoch ) = toc( start_epoch );

        %============================================
        % Collect cost function metrics, Plot results
        %============================================
        
        if ( mod( sol.it.epoch, sol.it.metrics_and_plotting ) == 0 ) || ( sol.it.epoch == 1 )

            sol.sample.T  = gather( sol.GPU.TF );
            sol.probe.phi = gather( sol.GPU.phi );
            
            [ sol ] = ptycho2DTPA_collectmetrics( sol, expt );
            
            ptycho2DTPA_plotresults( sol, expt );

        end

        %===============================
        % Feedback for iteration counter
        %===============================
        
        S2 = num2str( [ kk, N_epochs, sol.it.epoch, sol.timings.epoch( sol.it.epoch ) ], 'Local Epoch = %d / %d, Total Epochs = %d, t_epoch = %e' );

        fprintf( [ S2, '\n'])
        
        %========================
        % Epoch iteration counter
        %========================
                
        sol.it.epoch = sol.it.epoch + 1;  

    end
    
    sol.sample.T  = gather( sol.GPU.TF );
    sol.probe.phi = gather( sol.GPU.phi );
    
    sol = rmfield( sol, 'sample_sposview_indices' );
%     sol = rmfield( sol, 'psi' );
    sol = rmfield( sol, 'GPU' );

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    
    %========
    
    GPU.rPIE_alpha = gpuArray( sol.rPIE_alpha );
%     
%     GPU.RAAR_beta = gpuArray( sol.RAAR_beta );
% 
%     GPU.psi = gpuArray( sol.psi );  

    %========
    
    GPU.rs    = gpuArray( sol.spos.rs );
    GPU.Nspos = gpuArray( sol.spos.N );
    
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    GPU.sz      = gpuArray( sol.sz.sz );
    GPU.rc      = gpuArray( sol.sz.rc );
    GPU.sqrt_rc = gpuArray( sol.sz.sqrt_rc );

    GPU.samsz = gpuArray( sol.sample.sz.sz );
    GPU.samrc = gpuArray( sol.sample.sz.rc );

    %========

    GPU.ind        = uint32( sol.sample_sposview_indices ); 
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
    GPU.ind        = gpuArray( GPU.ind ); 
    GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %========
    
    GPU.meas_D = reshape( expt.meas.D, [ GPU.sz, 1, GPU.Nspos ] );

    GPU.meas_D    = gpuArray( GPU.meas_D );
    GPU.meas_Deq0 = gpuArray( GPU.meas_D == 0 );
    GPU.measLPF   = gpuArray( sol.measLPF );
    
    %========

    GPU.phi           = gpuArray( sol.probe.phi );
    GPU.fro2TOT       = gpuArray( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = gpuArray( sol.probe.scpm.occ );
    GPU.probe_support = gpuArray( sol.probe.support );
    GPU.scpmmax       = gpuArray( sol.probe.scpm.max );
    
    %========
    
    GPU.swparams.blurx     = gpuArray( sol.swparams.blurx ); 
    GPU.swparams.blury     = gpuArray( sol.swparams.blury );
    GPU.swparams.sparselvl = gpuArray( sol.swparams.sparselvl ); 
    
    %========
    
    GPU.TF         = gpuArray( sol.sample.T );
    GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);
    
    GPU.vs_r = gpuArray( sol.sample.vs.r );
    GPU.vs_c = gpuArray( sol.sample.vs.c );
    
    %========

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

    %========

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )


    %==========================================
    % Move relevant information onto GPU memory
    
    %========
    
    GPU.RAAR_beta = single( sol.RAAR.beta );
    
    GPU.psi = single( sol.psi );                                                       

    %========
    
    GPU.rs    = single( sol.spos.rs );
    GPU.Nspos = single( sol.spos.N );
    
    GPU.Nscpm = single( sol.probe.scpm.N );

    GPU.sz      = single( sol.sz.sz );
    GPU.rc      = single( sol.sz.rc );
    GPU.sqrt_rc = single( sol.sz.sqrt_rc );

    GPU.samsz = single( sol.sample.sz.sz );
    GPU.samrc = single( sol.sample.sz.rc );

    %========

    GPU.ind        = uint32( sol.sample_sposview_indices ); 
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
%     GPU.ind        = gpuArray( GPU.ind ); 
%     GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %========
    
    GPU.meas_D = reshape( expt.meas.D, [ GPU.sz, 1, GPU.Nspos ] );

    GPU.meas_D    = single( GPU.meas_D );
    GPU.meas_Deq0 = single( GPU.meas_D == 0 );
    GPU.measLPF   = single( sol.measLPF );
    
    %========

    GPU.phi           = single( sol.probe.phi );
    GPU.fro2TOT       = single( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = single( sol.probe.scpm.occ );
    GPU.probe_support = single( sol.probe.support );
    GPU.scpmmax       = single( sol.probe.scpm.max );
    
    %========
    
    GPU.swparams.blurx     = single( sol.swparams.blurx ); 
    GPU.swparams.blury     = single( sol.swparams.blury );
    GPU.swparams.sparselvl = single( sol.swparams.sparselvl ); 
    
    %========
    
    GPU.TF         = single( sol.sample.T );
    GPU.abs_TF_lim = single( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = single( [ sol.sample.phsL, sol.sample.phsH ]);
    
    GPU.vs_r = single( sol.sample.vs.r );
    GPU.vs_c = single( sol.sample.vs.c );
    
    %========

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

    %========

end

%====================================================================================================================================================
