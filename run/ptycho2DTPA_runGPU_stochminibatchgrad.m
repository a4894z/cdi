function [ sol, expt ] = ptycho2DTPA_runGPU_stochminibatchgrad( sol, expt, N_epochs )

    %==========================================
    % get name of current function being called
    %==========================================
    
    st = dbstack;
    namestr = st.name;

    %================================================================
    % Send to GPU parameters that will remain constant for all epochs
    %================================================================

    sol.GPU.rPIE_alpha_T   = gpuArray( sol.rPIE_alpha_T );
    sol.GPU.rPIE_alpha_phi = gpuArray( sol.rPIE_alpha_phi );

    %========

    sol.GPU.sample_sparse.threshname          = sol.sample.sparse.threshname;
    sol.GPU.sample_sparse.threshtype          = sol.sample.sparse.threshtype;
    sol.GPU.sample_sparse.lvl                 = gpuArray( sol.sample.sparse.lvl );
    sol.GPU.sample_sparse.support             = gpuArray( sol.sample.sparse.support ); 
    sol.GPU.sample_sparse.s2DFDxy.fft_scaling = gpuArray( sol.sample.sparse.s2DFDxy.fft_scaling );
    sol.GPU.sample_sparse.s2DFDxy.qLPF        = gpuArray( sol.sample.sparse.s2DFDxy.qLPF );
    sol.GPU.sample_sparse.s2DFDxy.Dx_fft      = gpuArray( sol.sample.sparse.s2DFDxy.Dx_fft );
    sol.GPU.sample_sparse.s2DFDxy.Dy_fft      = gpuArray( sol.sample.sparse.s2DFDxy.Dy_fft );
    sol.GPU.sample_sparse.s2DFDxy.Dx_fft_conj = gpuArray( sol.sample.sparse.s2DFDxy.Dx_fft_conj );
    sol.GPU.sample_sparse.s2DFDxy.Dy_fft_conj = gpuArray( sol.sample.sparse.s2DFDxy.Dy_fft_conj );
    sol.GPU.sample_sparse.s2DFDxy.sqrt_rc     = gpuArray( sol.sample.sparse.s2DFDxy.sqrt_rc );    
        
    %========

    sol.GPU.probe_sparse.threshtype = sol.probe.sparse.threshtype;
    sol.GPU.probe_sparse.threshname = sol.probe.sparse.threshname;
    
%     sol.GPU.probe_sparse.pct_x = gpuArray( sol.probe.sparse.pct_x );
%     sol.GPU.probe_sparse.pct_y = gpuArray( sol.probe.sparse.pct_y );
    sol.GPU.probe_sparse.lvl_x = gpuArray( sol.probe.sparse.lvl_x );
    sol.GPU.probe_sparse.lvl_y = gpuArray( sol.probe.sparse.lvl_y );
    
    sol.GPU.probe_sparse.support = gpuArray( sol.probe.sparse.support );
    
    sol.GPU.probe_sparse.s2DFDxy.fft_scaling = gpuArray( sol.probe.sparse.s2DFDxy.fft_scaling );
    sol.GPU.probe_sparse.s2DFDxy.qLPF        = gpuArray( sol.probe.sparse.s2DFDxy.qLPF );
    sol.GPU.probe_sparse.s2DFDxy.Dx_fft      = gpuArray( sol.probe.sparse.s2DFDxy.Dx_fft );
    sol.GPU.probe_sparse.s2DFDxy.Dy_fft      = gpuArray( sol.probe.sparse.s2DFDxy.Dy_fft );
    sol.GPU.probe_sparse.s2DFDxy.Dx_fft_conj = gpuArray( sol.probe.sparse.s2DFDxy.Dx_fft_conj );
    sol.GPU.probe_sparse.s2DFDxy.Dy_fft_conj = gpuArray( sol.probe.sparse.s2DFDxy.Dy_fft_conj );
    sol.GPU.probe_sparse.s2DFDxy.sqrt_rc     = gpuArray( sol.probe.sparse.s2DFDxy.sqrt_rc );

    %=================
    % probe shrinkwrap
    %=================
    
    sol.GPU.swparams.blurx     = gpuArray( sol.swparams.blurx ); 
    sol.GPU.swparams.blury     = gpuArray( sol.swparams.blury );
    sol.GPU.swparams.sparselvl = gpuArray( sol.swparams.sparselvl ); 
    
    %========    
     
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
    
    sol.GPU.TFvec = gpuArray( sol.sample.T( : ));
    
%     sol.GPU.psi           = gpuArray( sol.psi );
    
    sol.GPU.phi                = gpuArray( sol.probe.phi );
    sol.GPU.fro2TOT            = gpuArray( sol.probe.scpm.fro2TOT );
    sol.GPU.scpmocc            = gpuArray( sol.probe.scpm.occ );
    sol.GPU.probe_support      = gpuArray( sol.probe.support );
    sol.GPU.scpmmax            = gpuArray( sol.probe.scpm.max );
    sol.GPU.Nscpm              = gpuArray( sol.probe.scpm.N );
    sol.GPU.probe_gaussian_lpf = gpuArray( sol.probe.gaussian_lpf );
     
    %=================================================
    % step length parameters for Poisson cost function 
    %=================================================
    
    if strcmp( sol.exwv_noisemodel, 'poisson' )
        
        sol.GPU.Nalpha = gpuArray( single( 6 ));

        sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ]; ...
                                 [ 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 2.50, 2.50 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ]; ...
    %                              [ 6.00, 6.00, 4.00, 3.00, 3.00, 1.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ]; ...
    %                              [ 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00, 0.00, 0.00 ]; ...
    %                              [ 2.00, 2.00, 2.00, 2.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00, 0.00 ]; ...
    %                              [ 1e-6, 2.00, 2.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00, 0.00 ]; ...
    %                              [ 1.00, 1.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00, 0.00 ]; ...
    %                              [ 1.00, 1.00 ] ];

    %     sol.GPU.alpha_minmax = [ [ 0.00 ]; ...
    %                              [ 1.00 ] ];

        %========

        sol.GPU.alpha_minmax = gpuArray( single( sol.GPU.alpha_minmax ));

        for pp = 1 : sol.GPU.Nscpm

            sol.GPU.alpha_test( pp, : ) = linspace( sol.GPU.alpha_minmax( 1, pp ), sol.GPU.alpha_minmax( 2, pp ), sol.GPU.Nalpha );

        end

        sol.GPU.alpha_test = gpuArray( reshape( sol.GPU.alpha_test, [ 1, 1, sol.GPU.Nscpm, sol.GPU.Nalpha ] ));
        sol.GPU.alpha_prev = gpuArray( zeros( [ sol.spos.N, sol.GPU.Nscpm ], 'single' ) );
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        texwv   = 0;
        tsample = 0;
        tSCPM   = 0;
        
        meas_gauss_intensity = 0;
        meas_gauss_magnitude = 0;
        meas_poiss           = 0;
        
        grad_meas_gauss_intensity = 0;
        grad_meas_gauss_magnitude = 0;
        grad_meas_poiss           = 0;
        
        collect_metrics = (( mod( sol.it.epoch, sol.it.collect_metrics ) == 0 ) || ( sol.it.epoch == 1 ));
        
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exitwave Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            start_exwv = tic;

            %========

            if strcmp( sol.exwv_noisemodel, 'gaussian' )
            
                [ sol.GPU.psi, meas_metrics ] = exitwave_update_2DTPA_gaussian( sol.GPU.phi,       ...            
                                                                                sol.GPU.TFvec,     ...
                                                                                sol.GPU.ind,       ...
                                                                                sol.GPU.sz,        ...
                                                                                sol.GPU.Nspos,     ...
                                                                                sol.GPU.sqrt_rc,   ...
                                                                                sol.GPU.meas_D,    ...
                                                                                sol.GPU.meas_Deq0, ...
                                                                                sol.GPU.measLPF,   ...
                                                                                collect_metrics );
                                                                                                 
            elseif strcmp( sol.exwv_noisemodel, 'poisson' )
            
                if ( mod( sol.it.epoch, sol.it.exwv_poisson_exactsearch ) == 0 ) || ( sol.it.epoch == 1 )

                    [ sol.GPU.psi, sol.GPU.alpha_opt, meas_metrics ] = exitwave_update_2DTPA_poisson( sol.GPU.phi,        ...
                                                                                                      sol.GPU.TFvec,      ...
                                                                                                      sol.GPU.ind,        ...
                                                                                                      sol.GPU.sz,         ...
                                                                                                      sol.GPU.Nspos,      ...
                                                                                                      sol.GPU.Nscpm,      ...
                                                                                                      sol.GPU.alpha_test, ...
                                                                                                      [],                 ...
                                                                                                      sol.GPU.rc,         ... 
                                                                                                      sol.GPU.sqrt_rc,    ...
                                                                                                      sol.GPU.I_m,        ...
                                                                                                      sol.GPU.meas_Deq0,  ...
                                                                                                      sol.GPU.measLPF,    ...
                                                                                                      collect_metrics );
                                                                                                  
                    %======================================================================
                    % keep track of step length for the scan positions we just computed for
                    %======================================================================
                    
                    sol.GPU.alpha_prev( sol.GPU.batch_indx, : ) = squeeze( sol.GPU.alpha_opt );

                else

                    alpha_prev = sol.GPU.alpha_prev( sol.GPU.batch_indx, : );
                    alpha_prev = reshape( alpha_prev, [ 1, size( alpha_prev ) ] );

                    [ sol.GPU.psi, sol.GPU.alpha_opt, meas_metrics ] = exitwave_update_2DTPA_poisson( sol.GPU.phi,       ...
                                                                                                      sol.GPU.TFvec,     ...
                                                                                                      sol.GPU.ind,       ...
                                                                                                      sol.GPU.sz,        ...
                                                                                                      sol.GPU.Nspos,     ...
                                                                                                      sol.GPU.Nscpm,     ...
                                                                                                      [],                ...
                                                                                                      alpha_prev,        ...
                                                                                                      sol.GPU.rc,        ... 
                                                                                                      sol.GPU.sqrt_rc,   ...
                                                                                                      sol.GPU.I_m,       ...
                                                                                                      sol.GPU.meas_Deq0, ...
                                                                                                      sol.GPU.measLPF,   ...
                                                                                                      collect_metrics );

                end
            
            else
        
                error( 'Invalid Exitwave Update Noise Model: Use Only ''gaussian'' or ''poisson'' (Case Sensitive)' )
            
            end

            %=======
            
            if collect_metrics
                
                meas_gauss_intensity = meas_gauss_intensity + meas_metrics.gauss_intensity; 
                meas_gauss_magnitude = meas_gauss_magnitude + meas_metrics.gauss_magnitude; 
                meas_poiss           = meas_poiss           + meas_metrics.poiss; 
                
                grad_meas_gauss_intensity = grad_meas_gauss_intensity + meas_metrics.grad_gauss_intensity; 
                grad_meas_gauss_magnitude = grad_meas_gauss_magnitude + meas_metrics.grad_gauss_magnitude; 
                grad_meas_poiss           = grad_meas_poiss           + meas_metrics.grad_poiss; 
                
            end
   
            texwv = texwv + toc( start_exwv );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sample Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ( mod( sol.it.epoch, sol.it.sample_update ) == 0 ) && ( sol.it.epoch >= sol.it.sample_start )

                start_sample = tic;

                %========================================================================
                % keep a copy of the current sample for use in probe update testing below
                %========================================================================

                sol.GPU.TFvec_old = sol.GPU.TFvec;

                %========================================
                % Vectorized minibatch rPIE sample update
                %========================================
            
                [ sol.GPU.TFvec ] = rPIEupdate_batch_2DTPA_sample( sol.GPU.psi,        ...
                                                                   sol.GPU.TFvec,      ...
                                                                   sol.GPU.phi,        ...
                                                                   sol.GPU.ind_offset, ...
                                                                   sol.GPU.rc,         ...
                                                                   sol.GPU.Nspos,      ...
                                                                   sol.GPU.rPIE_alpha_T );           
                                                              
                %================
                % Sample supports
                %================

                if mod( sol.it.epoch, sol.it.sample_sparsity ) == 0

                    %==========================
                    % Shrinkwrap sample support
                    %==========================

        %             TF  = reshape( sol.GPU.TFvec, sol.GPU.samsz );
        %             TF  = lpf_gauss( TF, 0.70 * sol.GPU.samsz );
        %             TF  = sparseFDxy_update_2Dsample( TF, sol.GPU.sample_sparse );      
        %             sol.GPU.TFvec = TF( : );   
        %             clear( 'TF' )

                    %===================================================================
                    % Sample analysis sparsity on reduced dimensionality representations
                    %===================================================================

                    [ sol.GPU.TFvec ] = sparseFDxy_update_2Dsample( reshape( sol.GPU.TFvec, sol.GPU.samsz ), sol.GPU.sample_sparse );

                    sol.GPU.TFvec = sol.GPU.TFvec( : );

                end

                %========================================
                % Sample mag/phase inequality constraints 
                %========================================

                if mod( sol.it.epoch, sol.it.sample_mag_ineq ) == 0

                    sol.GPU.TFvec = modulus_limits_project( sol.GPU.TFvec, sol.GPU.abs_TF_lim );
        %             sol.GPU.TFvec = modulus_limits_scale( sol.GPU.TFvec, sol.GPU.abs_TF_lim );

                end
                
                if mod( sol.it.epoch, sol.it.sample_phs_ineq ) == 0
                    
                    sol.GPU.TFvec = phase_limits_project( sol.GPU.TFvec, sol.GPU.phs_TF_lim );
%                     sol.GPU.TFvec = phase_limits_scale( sol.GPU.TFvec, sol.GPU.phs_TF_lim );

                end
                
                %======================================================
                % Sample mag and phase correlation ( weighted average )
                %======================================================

%                     if 0 %mod( sol.it.epoch, 50 ) == 0
% 
%                         abs_TF = abs( sol.GPU.TFvec );
%                         abs_TF = abs_TF - min( abs_TF( : ));
%                         abs_TF = abs_TF / max( abs_TF( : ));
% 
%                         phs_TF = angle( sol.GPU.TFvec );
%                         phs_TF = phs_TF - min( phs_TF( : ));
%                         phs_TF = phs_TF / max( phs_TF( : ));
% 
%             %             as = 0.75;
%             %             ps = 0.9;
%             %             sol.GPU.TFvec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
% 
%                         as = 0.5;
%                         sol.GPU.TFvec = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.TFvec ));
% 
%                     end  
                
                tsample = tsample + toc( start_sample );

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Probe Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch >= sol.it.probe_start )

                start_probe = tic;

                %=======================================
                % Vectorized minibatch rPIE probe update
                %=======================================
  
                T_view = reshape( sol.GPU.TFvec_old( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
     
                z   = sum( abs( T_view ) .^ 2, 4 );
                w_T = sol.GPU.rPIE_alpha_phi * ( max( z( : )) - z );
                
                sol.GPU.phi = ( sum( conj( T_view ) .* sol.GPU.psi, 4 ) + w_T .* sol.GPU.phi ) ./ ( z + w_T );

                %======================
                % Blur/smooth the SCPMs
                %======================
                
                if ( mod( sol.it.epoch, sol.it.probe_smoothing ) == 0 ) 
                    
                    tmp0 = fftshift( fftshift( sol.GPU.phi, 1 ), 2 );

                    tmp0 = fft( fft( tmp0, [], 1 ), [], 2 )   / sol.GPU.sqrt_rc;
                    tmp0 = tmp0 .* sol.GPU.probe_gaussian_lpf;
                    tmp0 = ifft( ifft( tmp0, [], 1 ), [], 2 ) * sol.GPU.sqrt_rc;

                    sol.GPU.phi = fftshift( fftshift( tmp0, 1 ), 2 );

                end
                
                %==========================================================================
                % Center the probe to middle of its array; shift sample as well accordingly
                %==========================================================================

                if ( mod( sol.it.epoch, sol.it.probe_centering ) == 0 ) 

%                     [~, Ic ] = max( sum( abs( sol.GPU.phi( :, :, end )), 1 ));
%                     [~, Ir ] = max( sum( abs( sol.GPU.phi( :, :, end )), 2 ));
%                     shift_px = gather( -1 * round( [ Ir, Ic ] - 0.5 * sol.GPU.sz - 1 ));

                    [ com ] = centerofmass( abs( sol.GPU.phi( :, :, end )));
%                     [ com ] = centerofmass( sum( abs( sol.GPU.phi .^ 2 ), 3 ));
                    shift_px = gather( -1 * round( com - 0.5 * sol.GPU.sz - 1 ) );
                    
                    if sum( shift_px ) ~= 0

                        for pp = 1 : sol.GPU.Nscpm

                            sol.GPU.phi( :, :, pp ) = nocircshift2D( sol.GPU.phi( :, :, pp ), shift_px );

                        end
                        
                        sol.GPU.TFvec = nocircshift2D( reshape( sol.GPU.TFvec, sol.GPU.samsz ), shift_px );
                        sol.GPU.TFvec = sol.GPU.TFvec( : );

                        sol.GPU.TFvec( sol.GPU.TFvec == 0 ) = 1;
                        
                    end
                    
                end

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

%                     tmp0          = sum( abs( sol.GPU.phi ) .^ 2, 3 );
%                     [ ~, supp ]   = shrinkwrap( tmp0, sol.GPU.swparams ); 
%                     sol.GPU.phi   = sol.GPU.phi .* supp;

                    %==================================
                    % Support using dominant probe mode
                    %==================================

                    [ sol.GPU.phi( :, :, end ), supp ] = shrinkwrap( sol.GPU.phi( :, :, end ), sol.GPU.swparams ); 
                    sol.GPU.phi = sol.GPU.phi .* supp;

                end

                %===================
                % Max abs constraint
                %===================

                if ( mod( sol.it.epoch, sol.it.probe_maxvals ) == 0 ) && ~isempty( sol.GPU.scpmmax )
                    
                    tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
%                     tmp2 = sol.GPU.scpmmax;
                    tmp0 = ( abs( sol.GPU.phi ) > tmp2 );
                    tmp1 = not( tmp0 );
                        
                    sol.GPU.phi = sol.GPU.phi .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.phi )) .* tmp0;
                    
                    clear( 'tmp0', 'tmp1' )
                    
                end

                %===========================================================
                % Probe Scaling ( # Photons ) and SCPM occupancy constraints
                %===========================================================

                if ( mod( sol.it.epoch, sol.it.probe_scaling ) == 0 ) 
                    
                    [ sol.GPU.phi, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.phi, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

                end

                %==================================
                % TESTING: probe phase ramp removal
                %==================================
                
                if ( mod( sol.it.epoch, sol.it.probe_phase_sparse ) == 0 ) 
                    
                    for pp = 1 : sol.GPU.Nscpm
                        
                        phase_phi = angle( sol.GPU.phi( :, :, pp ) );

                        [ edges_phase_phi ] = edgedetect_FDxy( phase_phi, sol.GPU.probe_sparse.s2DFDxy );

                        abs_edges_phase_phi.y = abs( edges_phase_phi.y );
                        abs_edges_phase_phi.x = abs( edges_phase_phi.x );

                        lamr = find_thresh_from_sparsitylevel( abs_edges_phase_phi.y, sol.GPU.probe_sparse.lvl_y );
                        lamc = find_thresh_from_sparsitylevel( abs_edges_phase_phi.x, sol.GPU.probe_sparse.lvl_x );

                        [ edges_phase_phi.y, Wr ] = soft_shrinkage( edges_phase_phi.y, abs_edges_phase_phi.y, lamr );
                        [ edges_phase_phi.x, Wc ] = soft_shrinkage( edges_phase_phi.x, abs_edges_phase_phi.x, lamc );     

                        [ phase_phi ] = iedgedetect_FDxy( edges_phase_phi, sol.GPU.probe_sparse.s2DFDxy );

                        sol.GPU.phi( :, :, pp ) = abs( sol.GPU.phi( :, :, pp ) ) .* exp( 1i * real( phase_phi ));     
                    
                    end
                    
                    
%                     % TESTING: REAL VALUED PROBE
%                     sol.GPU.phi = abs( sol.GPU.phi );

                end
                

                %========
                
                tSCPM = tSCPM + toc( start_probe );
   
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Batch updates complete, orthogonalize final SCPMs for this epoch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ( mod( sol.it.epoch, sol.it.probe_update ) == 0 ) && ( sol.it.epoch >= sol.it.probe_start )
        
            start_probe = tic;
            
            %========================================================
            % Orthogonalize the SCPMs (spatial coherence probe modes)
            %========================================================

            if ( mod( sol.it.epoch, sol.it.probe_orthog ) == 0 ) || ( sol.it.epoch == 1 )

                [ sol.GPU.phi ] = orthog_modes_eigendecomp( sol.GPU.phi ); 
                
                % get new occupancies resulting from orthogonalization
                [ scpm ] = compute_scpm_photonocc( sol.GPU.phi );

                sol.GPU.fro2TOT = scpm.fro2TOT;
                sol.GPU.scpmocc = scpm.occ;

            end
                
            tSCPM = tSCPM + toc( start_probe );
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics and Misc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %===============
        % Timing Metrics
        %===============

        sol.timings.exwv_update(   sol.it.epoch ) = texwv;        
        sol.timings.sample_update( sol.it.epoch ) = tsample;
        sol.timings.probe_update(  sol.it.epoch ) = tSCPM;
        sol.timings.epoch(         sol.it.epoch ) = toc( start_epoch );
   
        %================================================
        % Measurement cost function metrics, plot results
        %================================================
        
        if collect_metrics 
            
            sol.metrics.meas_gauss_intensity( sol.it.metr ) = gather( meas_gauss_intensity ) / sol.spos.N;
            sol.metrics.meas_gauss_magnitude( sol.it.metr ) = gather( meas_gauss_magnitude ) / sol.spos.N;
            sol.metrics.meas_poiss(           sol.it.metr ) = gather( meas_poiss )           / sol.spos.N;
            
            sol.metrics.grad_meas_gauss_intensity( sol.it.metr ) = gather( grad_meas_gauss_intensity ) / ( sol.probe.scpm.N * sol.spos.N );
            sol.metrics.grad_meas_gauss_magnitude( sol.it.metr ) = gather( grad_meas_gauss_magnitude ) / ( sol.probe.scpm.N * sol.spos.N );
            sol.metrics.grad_meas_poiss(           sol.it.metr ) = gather( grad_meas_poiss )           / ( sol.probe.scpm.N * sol.spos.N );
            
            sol.it.mtot( sol.it.metr ) = sol.it.epoch;   
            sol.it.metr                = sol.it.metr + 1;   
            
        end
        
        %========
        
        if ( mod( sol.it.epoch, sol.it.mkimg_meas_metric ) == 0 )
            
            ptycho2DTPA_mkimg_meas_metric( sol, expt ); 
        
        end
        
        %========
        
        if ( mod( sol.it.epoch, sol.it.mkimg_sample_SCPM ) == 0 )
            
            sol.sample.T       = gather( reshape( sol.GPU.TFvec, sol.sample.sz.sz ));
            sol.probe.phi      = gather( sol.GPU.phi );
%             sol.probe.scpm.occ = gather( sol.GPU.scpmocc );
            
            ptycho2DTPA_mkimg_sample_SCPM( sol, expt ); 
            
        end

        %===============================
        % Feedback for iteration counter
        %===============================
        
        S2 = num2str( [ kk, ...
                        N_epochs, ...
                        sol.it.epoch, ...
                        sol.spos.rand_spos_subset_pct, ...
                        sol.timings.epoch( sol.it.epoch ) ], ...
                        'Local Epoch = %d / %d, Total Epochs = %d, minibatch pct = %0.4f, t_epoch = %e' );

        fprintf( [ S2, '\n'])
        
        if mod( kk, 10 ) == 0
            
            fprintf( [ '\n', namestr, '\n', pwd, '\n', 'Data mat name = ', expt.paths.rsdata, '\n\n' ])

        end
        
        %========================
        % Epoch iteration counter
        %========================
                
        sol.it.epoch = sol.it.epoch + 1;  
                
    end
    
    sol.sample.T       = gather( reshape( sol.GPU.TFvec, [ sol.sample.sz.sz ]));
    sol.probe.phi      = gather( sol.GPU.phi );
    sol.probe.scpm.occ = gather( sol.GPU.scpmocc );
    
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
 
    if strcmp( sol.exwv_noisemodel, 'gaussian' )
        
        %==========================================================================================
        % measurement amplitudes for when we use Gaussian cost function w/ constant (ignored) stdev
        %==========================================================================================

        GPU.meas_D = gpuArray( reshape( expt.meas.D( :, :, GPU.batch_indx ), [ GPU.sz, 1, GPU.Nspos ] ));
        GPU.meas_Deq0 = ( GPU.meas_D == 0 );

    elseif strcmp( sol.exwv_noisemodel, 'poisson' )
        
        %==============================================================
        % measurement intensities for when we use Poisson cost function
        %==============================================================

        GPU.I_m = gpuArray( reshape( expt.meas.D( :, :, GPU.batch_indx ) .^ 2, [ GPU.rc, GPU.Nspos ] ));
        GPU.meas_Deq0 = ( GPU.I_m == 0 );

    else
        
        error( 'Invalid Exitwave Update Noise Model: Use Only ''gaussian'' or ''poisson'' (Case Sensitive)' )
            
    end
    
end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )

    GPU = [];

end
