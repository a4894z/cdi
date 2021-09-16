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
        
        start_epoch = tic;
                
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

            for gg = 1 : N_repeat

                %====================================================================================================================================
                % Exitwave Update
                %================
                
                start_exwv = tic;
                
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
    
                sol.timings.exwv_update( sol.it.epoch ) = toc( start_exwv );
   
                %====================================================================================================================================
                % Sample Update
                %==============

                
                
                
                
                
                
                
                
                
                
%                 %======================================================================================
%                 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 % EXPERIMENTAL: keep a copy of the current sample for use in probe update testing below
%                 
%                 sol.GPU.T_vec_old = sol.GPU.T_vec;
%                 
%                 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 %======================================================================================
                
                
                
                
                
                
                
                
                
                start_sample = tic;
                
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
                
                sol.timings.sample_update( sol.it.epoch ) = toc( start_sample );
                
                
                
                
                
                
                
                
                
                
                
%                 %================================================================
%                 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 % EXPERIMENTAL: another exitwave update using just updated sample
%  
%                 sol.GPU.psi = ER_GPU_arrays_hadamard( sol.GPU.phi,       ...            % ERvec_2DTPA
%                                                       sol.GPU.T_vec,     ...
%                                                       sol.GPU.ind,       ...
%                                                       sol.GPU.sz,        ...
%                                                       sol.GPU.Nspos,     ...
%                                                       sol.GPU.sqrt_rc,   ...
%                                                       sol.GPU.meas_D,    ...
%                                                       sol.GPU.meas_Deq0, ...
%                                                       sol.GPU.measLPF );
%                                                   
%                 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 %================================================================
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                %====================================================================================================================================
                % Probe Update
                %=============
                
                start_probe = tic;
                
                if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch > sol.it.probe_start )

                    
                    
                    %==================================
                    % DM/ePIE hybrid style probe update  
                    %==================================

                    T_view = reshape( sol.GPU.T_vec( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);

                    abs2_TFview = abs( T_view ) .^ 2;

                    sol.GPU.phi = sol.GPU.phi + ( sum( conj( T_view ) .* sol.GPU.psi - sol.GPU.phi .* abs2_TFview, 4 ) ) ./ ( 1e-7 + sum( abs2_TFview, 4 ) );
                    

                    
                    
                    
                    
                    
                    
                    

%                     %======================================================
%                     %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                     % EXPERIMENTAL: use old sample estimate to update probe
% 
%                     T_view = reshape( sol.GPU.T_vec_old( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
% 
%                     abs2_TFview = abs( T_view ) .^ 2;
% 
%                     sol.GPU.phi = sol.GPU.phi + ( sum( conj( T_view ) .* sol.GPU.psi - sol.GPU.phi .* abs2_TFview, 4 ) ) ./ ( 1e-7 + sum( abs2_TFview, 4 ) );
% 
%                     %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                     %======================================================

                
                
                
                
                
                
                
                
                

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

                    %========================================================
                    % Orthogonalize the SCPMs (spatial coherence probe modes)
                    %========================================================

                    if ( mod( sol.it.epoch, sol.it.probe_orthog ) == 0 ) || ( sol.it.epoch == 1 )

                        [ sol.GPU.phi ] = orthog_modes_eigendecomp( sol.GPU.phi ); 

                    end
        
                end
                
                sol.timings.probe_update( sol.it.epoch ) = toc( start_probe );
                
            end

        end

        %==============
        % Epoch timing
        %=============
        
        sol.timings.epoch( sol.it.epoch ) = toc( start_epoch );
   
        %============================================
        % Collect cost function metrics, Plot results
        %============================================
        
        if ( mod( sol.it.epoch, sol.it.metrics_and_plotting ) == 0 ) || ( sol.it.epoch == 1 )

            sol.sample.T  = gather( reshape( sol.GPU.T_vec, [ sol.sample.sz.sz ]));
            sol.probe.phi = gather( sol.GPU.phi );
            
            [ sol ] = ptycho2DTPA_collectmetrics( sol, expt );
            
            ptycho2DTPA_plotresults( sol, expt );

        end
        
        %===============================
        % Feedback for iteration counter
        %===============================
        
        S2 = num2str( [ kk, N_epochs, sol.it.epoch, sol.spos.rand_spos_subset_pct, sol.timings.epoch( sol.it.epoch ) ], 'Local Epoch = %d / %d, Total Epochs = %d, minibatch pct = %0.4f, t_epoch = %e' );

        fprintf( [ S2, '\n'])
        
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
