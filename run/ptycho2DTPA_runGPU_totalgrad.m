function [ sol, expt ] = ptycho2DTPA_runGPU_totalgrad( sol, expt, N_epochs )

    %==========================================
    % get name of current function being called
    %==========================================
    
    st = dbstack;
    namestr = st.name;
    
    %=======================================================
    % form initial exitwaves from current sample/probe modes
    %=======================================================
    
    sol.spos.batch_indxsubset = randperm( sol.spos.N );
%     sol.spos.batch_indxsubset = 1 : sol.spos.N;
    
    sol.spos.batch_rs = sol.spos.rs( sol.spos.batch_indxsubset, : );
    
    sol.sample_sposview_indices = get_indices_2Dframes( sol.spos.batch_rs, sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );

%     sol.psi = sol.probe.phi .* reshape( sol.sample.T( sol.sample_sposview_indices ), [ sol.sz.sz, 1, size( sol.sample_sposview_indices, 2 ) ]);

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
        
        %========
        
        collect_metrics = (( mod( sol.it.epoch, sol.it.collect_metrics ) == 0 ) || ( sol.it.epoch == 1 ));
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exitwave Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        start_exwv = tic;
        
        %========

        [ sol.GPU.psi, meas_L2metric ] = exitwave_update_2DTPA_gaussian( sol.GPU.phi,       ...                      
                                                                         sol.GPU.TFvec,     ...
                                                                         sol.GPU.ind,       ...
                                                                         sol.GPU.sz,        ...
                                                                         sol.GPU.Nspos,     ...
                                                                         sol.GPU.sqrt_rc,   ...
                                                                         sol.GPU.meas_D,    ...
                                                                         sol.GPU.meas_Deq0, ...
                                                                         sol.GPU.measLPF,   ...
                                                                         collect_metrics );

        %========
        
        sol.timings.exwv_update( sol.it.epoch ) = toc( start_exwv );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sample Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if mod( sol.it.epoch, sol.it.sample_update ) == 0
            
            start_sample = tic;

            %========================================================================       
            % keep a copy of the current sample for use in probe update testing below
            %========================================================================       
            
            sol.GPU.TFvec_old = sol.GPU.TFvec;
            
            %========

            [ sol.GPU.TFvec ] = rPIEupdate_batch_2DTPA_sample( sol.GPU.psi,          ...
                                                               sol.GPU.TFvec,        ...
                                                               sol.GPU.phi,          ...
                                                               sol.GPU.ind_offset,   ...
                                                               sol.GPU.rc,           ...
                                                               sol.GPU.Nspos,        ...
                                                               sol.GPU.rPIE_alpha_T );                               

            %====================
            % Sparse Sample Edges
            %====================

            if mod( sol.it.epoch, sol.it.sample_sparsity ) == 0

                %===========================
                % Sparsity using TV variants
                %===========================

                [ sol.GPU.TFvec ] = sparseFDxy_update_2Dsample( reshape( sol.GPU.TFvec, sol.GPU.samsz ), sol.GPU.sparseGPU );

                sol.GPU.TFvec = sol.GPU.TFvec( : );

            end

            %======================================================
            % Sample inequality constraints (projection operations)
            %======================================================

%             if mod( sol.it.epoch, sol.it.sample_mag_phs_ineq ) == 0
% 
%                 sol.GPU.TFvec = modulus_limits_project( sol.GPU.TFvec, sol.GPU.abs_TF_lim );
% 
%     %             sol.GPU.TFvec = modulus_limits_project( sol.GPU.TFvec, [ 0, 1 ] );
%     %             sol.GPU.TFvec = modulus_limits_scale( sol.GPU.TFvec, sol.GPU.abs_TF_lim );
% 
% %                     sol.GPU.TFvec = phase_limits_project( sol.GPU.TFvec, sol.GPU.phs_TF_lim );
% %                     sol.GPU.TFvec = phase_limits_scale( sol.GPU.TFvec, sol.GPU.phs_TF_lim );
% 
%             end
            
            
            if mod( sol.it.epoch, sol.it.sample_mag_ineq ) == 0

                sol.GPU.TFvec = modulus_limits_project( sol.GPU.TFvec, sol.GPU.abs_TF_lim );
%                 sol.GPU.TFvec = modulus_limits_scale( sol.GPU.TFvec, sol.GPU.abs_TF_lim );

            end

            if mod( sol.it.epoch, sol.it.sample_phs_ineq ) == 0

                sol.GPU.TFvec = phase_limits_project( sol.GPU.TFvec, sol.GPU.phs_TF_lim );
%                 sol.GPU.TFvec = phase_limits_scale( sol.GPU.TFvec, sol.GPU.phs_TF_lim );

            end

            %======================================================
            % Sample mag and phase correlation ( weighted average )
            %======================================================

%             if 0 %mod( sol.it.epoch, 50 ) == 0
% 
%                 abs_TF = abs( sol.GPU.TFvec );
%                 abs_TF = abs_TF - min( abs_TF( : ));
%                 abs_TF = abs_TF / max( abs_TF( : ));
% 
%                 phs_TF = angle( sol.GPU.TFvec );
%                 phs_TF = phs_TF - min( phs_TF( : ));
%                 phs_TF = phs_TF / max( phs_TF( : ));
% 
%     %             as = 0.75;
%     %             ps = 0.9;
%     %             sol.GPU.TFvec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
% 
%                 as = 0.5;
%                 sol.GPU.TFvec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.TFvec ));
% 
%             end  

            sol.timings.sample_update( sol.it.epoch ) = toc( start_sample );
            
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Probe Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch > sol.it.probe_start )

            start_probe = tic;

            %=============================
            % Vectorized rPIE probe update
            %=============================
 
            T_view = reshape( sol.GPU.TFvec_old( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);

            z   = sum( abs( T_view ) .^ 2, 4 );
            w_T = sol.GPU.rPIE_alpha_phi * ( max( z( : )) - z );

            sol.GPU.phi = ( sum( conj( T_view ) .* sol.GPU.psi, 4 ) + w_T .* sol.GPU.phi ) ./ ( z + w_T );
 
            %==============
            % Probe Support
            %==============

            if ( mod( sol.it.epoch, sol.it.probe_support ) == 0 ) 
                
                %=========================
                % Fixed predefined support
                %=========================

                sol.GPU.phi = sol.GPU.phi .* sol.GPU.probe_support; 

               %=============================================
                % Shrinkwrap support from probe mode intensity
                %=============================================

                tmp0        = sum( abs( sol.GPU.phi ) .^ 2, 3 );
                [ ~, supp ] = shrinkwrap( tmp0, sol.GPU.swparams ); 
                sol.GPU.phi = sol.GPU.phi .* supp;

                %==================================
                % Support using dominant probe mode
                %==================================

%                     [ sol.GPU.phi( :, :, 3 ), supp ] = shrinkwrap( sol.GPU.phi( :, :, 3 ), swparams ); 
%                     sol.GPU.phi = sol.GPU.phi .* supp;

            end

            %===================
            % Max abs constraint
            %===================

            if ( mod( sol.it.epoch, sol.it.probe_maxvals ) == 0 ) && ~isempty( sol.GPU.scpmmax )

%                     tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
                tmp2 = sol.GPU.scpmmax;
                tmp0 = ( abs( sol.GPU.phi ) > tmp2 );
                tmp1 = not( tmp0 );

                sol.GPU.phi = sol.GPU.phi .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.phi )) .* tmp0;

                clear( 'tmp0', 'tmp1' )

            end

%             %===================================
%             % Center SCPMs using probe intensity
%             %===================================
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
            
            %===========================================================
            % Probe scaling ( # Photons ) and SCPM occupancy constraints
            %===========================================================

            if ( mod( sol.it.epoch, sol.it.probe_scaling ) == 0 ) || ( sol.it.epoch == 1 )

                [ sol.GPU.phi, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.phi, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

            end

        
            sol.timings.probe_update( sol.it.epoch ) = toc( start_probe );

        end

        %=============
        % Epoch timing
        %=============
        
        sol.timings.epoch( sol.it.epoch ) = toc( start_epoch );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics and Misc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if collect_metrics 
            
            sol.metrics.meas_all( sol.it.metr ) = gather( meas_L2metric ) / sol.spos.N;
            sol.it.mtot( sol.it.metr )          = sol.it.epoch;   
            sol.it.metr                         = sol.it.metr + 1;   
            
        end

        if ( mod( sol.it.epoch, sol.it.mkimg_meas_metric ) == 0 ), ptycho2DTPA_mkimg_meas_metric( sol, expt ); end
        
        if ( mod( sol.it.epoch, sol.it.mkimg_sample_SCPM ) == 0 )
            
            sol.sample.T  = gather( reshape( sol.GPU.TFvec, sol.sample.sz.sz ));
            sol.probe.phi = gather( sol.GPU.phi );

            ptycho2DTPA_mkimg_sample_SCPM( sol, expt ); 
            
        end
        
%         if ( mod( sol.it.epoch, sol.it.collect_metrics ) == 0 ) || ( sol.it.epoch == 1 )
% 
%             sol.sample.T  = gather( reshape( sol.GPU.TFvec, sol.sample.sz.sz ));
%             sol.probe.phi = gather( sol.GPU.phi );
%             
%             [ sol ] = ptycho2DTPA_collectmetrics( sol, expt );     
%             
%             if ( mod( sol.it.epoch, sol.it.sample_probe_img ) == 0 )
%                 
%                 ptycho2DTPA_plotresults( sol, expt );
%             
%             end
% 
%         end
        
        %===============================
        % Feedback for iteration counter
        %===============================
        
        S2 = num2str( [ kk, N_epochs, sol.it.epoch, sol.timings.epoch( sol.it.epoch ) ], 'Local Epoch = %d / %d, Total Epochs = %d, t_epoch = %e' );

        fprintf( [ S2, '\n'])

        if mod( kk, 10 ) == 0
            
            fprintf( [ '\n', namestr, '\n', pwd, '\n', 'Data mat name = ', expt.paths.rsdata, '\n\n' ])

        end
        
        %========================
        % Epoch iteration counter
        %========================
                
        sol.it.epoch = sol.it.epoch + 1;  

    end
    
    sol.sample.T  = gather( reshape( sol.GPU.TFvec, [ sol.sample.sz.sz ]));
    sol.probe.phi = gather( sol.GPU.phi );
    
    sol = rmfield( sol, 'sample_sposview_indices' );
%     sol = rmfield( sol, 'psi' );
    sol = rmfield( sol, 'GPU' );

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )
    
    GPU.device = gpuDevice;
    
    %==========================================
    % Move relevant information onto GPU memory
    %==========================================
    
%     GPU.RAAR_beta = gpuArray( single( sol.RAAR.beta ));
    
    GPU.rPIE_alpha_T   = gpuArray( sol.rPIE_alpha_T );
    GPU.rPIE_alpha_phi = gpuArray( sol.rPIE_alpha_phi );
    
    %========
    
    GPU.batch_N          = gpuArray( uint32( length( sol.spos.batch_indxsubset )));
    GPU.batch_indxsubset = gpuArray( uint32( sol.spos.batch_indxsubset ));
    GPU.batch_rs         = gpuArray( single( sol.spos.batch_rs ));
    
    GPU.Nspos = gpuArray( single( GPU.batch_N ));
    GPU.Nscpm = gpuArray( single( sol.probe.scpm.N ));

    %========
    
    GPU.sz      = gpuArray( single( sol.sz.sz ));
    GPU.rc      = gpuArray( single( sol.sz.rc ));
    GPU.sqrt_rc = gpuArray( single( sol.sz.sqrt_rc ));

    %========

    GPU.samsz = gpuArray( sol.sample.sz.sz );
    GPU.samrc = gpuArray( sol.sample.sz.rc );
    
    %========

    GPU.ind = uint32( sol.sample_sposview_indices ); 

%     GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
    GPU.ind        = gpuArray( GPU.ind ); 
    GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %========

    % meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = reshape( expt.meas.D( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
    GPU.meas_D = gpuArray( GPU.meas_D );
    
    GPU.meas_Deq0 = ( GPU.meas_D == 0 );

    GPU.measLPF = gpuArray( sol.measLPF );
    
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
    
    GPU.TFvec      = gpuArray( sol.sample.T( : ));
    GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);

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

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    %==========================================
    
    GPU.RAAR_beta = single( sol.RAAR.beta );
    
    %==============
    
    GPU.batch_N          = single( length( sol.spos.batch_indxsubset ) );
    GPU.batch_indxsubset = single( sol.spos.batch_indxsubset );
    GPU.batch_rs         = single( sol.spos.batch_rs );
    
    GPU.Nspos = single( GPU.batch_N );
    GPU.Nscpm = single( sol.probe.scpm.N );

    %==============
    
    GPU.sz      = single( sol.sz.sz );
    GPU.rc      = single( sol.sz.rc );
    GPU.sqrt_rc = single( sol.sz.sqrt_rc );

    %==============

    GPU.samsz = single( sol.sample.sz.sz );
    GPU.samrc = single( sol.sample.sz.rc );

    %==============
    
    GPU.psi = single( sol.psi );                                                       

    %==============

    GPU.ind = uint32( sol.sample_sposview_indices ); 

%     GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
%     GPU.ind        = gpuArray( GPU.ind ); 
%     GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %==============

    % meas_D = single( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = reshape( expt.meas.D( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
    GPU.meas_D = single( GPU.meas_D );
    
    GPU.meas_Deq0 = ( GPU.meas_D == 0 );

    GPU.measLPF = reshape( sol.measLPF, [ GPU.sz, 1, 1 ] );
    GPU.measLPF = single( sol.measLPF );
    
    %==============

    GPU.phi           = single( sol.probe.phi );
    GPU.fro2TOT       = single( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = single( sol.probe.scpm.occ );
    GPU.probe_support = single( sol.probe.support );
%     GPU.scpmmax       = single( sol.probe.scpm.max );
    
    %==============
    
    GPU.swparams.blurx     = single( sol.swparams.blurx ); 
    GPU.swparams.blury     = single( sol.swparams.blury );
    GPU.swparams.sparselvl = single( sol.swparams.sparselvl ); 
    
    %==============
    
    GPU.TFvec      = single( sol.sample.T( : ));
    GPU.abs_TF_lim = single( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = single( [ sol.sample.phsL, sol.sample.phsH ]);

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

%     % parallel.gpu.single.zeros(size(data), 'single');
% 
%     GPU.spos_conjP_exwv = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
%     GPU.spos_abs_P_abs2 = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));

    %====================================================
    % Ensure all relevant information is single precision
    %====================================================
    
    % ???

end
