function [ sol, expt ] = ptycho2DTPA_runGPUiterations_v3( sol, expt, N_datamovements, N_exwvup )

    %==========
    % Main Loop
    %==========
 
    for kk = 1 : N_datamovements

        %===================================================================
        % get new random subset of scan position indices for exitwave update 
        %===================================================================
        
        sol.spos.batch_indxsubset = randperm( length( sol.spos.indxsubset ), sol.spos.batch_N );
        sol.spos.batch_rs         = sol.spos.rs( sol.spos.batch_indxsubset, : );
        sol.spos.frameindx        = get_indices_2Dframes( sol.spos.batch_rs, sol.sample.sz.sz, sol.sample.vs );

        % figure; 
        % plot_2Dscan_positions( sol.spos.rs, [], sol.spos.batch_rs, [] )
        % %     plot_2Dscan_positions( expt.spos.yv_xh_frameL, [], expt.spos.yv_xh_frameL( 1 : 1 : end, : ), expt.spos.indxsubset( 1 : 1 : end ))
        % set( gca, 'xdir', 'reverse' )
        % set( gca, 'ydir', 'normal' )
        % grid on
        % %daspect([1 1 1])
        % xlabel('xh, lab frame'); ylabel('yv, lab frame');
        % %title('positions for scanning the probe on the sample')

        %=======================================================
        % form initial exitwaves from current sample/probe modes
        %=======================================================
        
        TF = sol.sample.TF( : );
        TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, sol.spos.batch_N ]);
        
        sol.phi = TF .* sol.probe.P;
        
        clear( 'TF' )

        %=========================================================
        % perform data movement from CPU main memory to GPU memory
        %=========================================================
        
        [ sol.GPU ] = CPUmem2GPUmem( sol, expt );

        %===============================================
        % perform exitwave, sample, probe updates on GPU
        %===============================================
        
%         exwv_diff = zeros( 1, N_exwvup, 'single' );

        for gg = 1 : N_exwvup
    
            %===============================
            % Feedback for iteration counter
            %===============================
            
            fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
            if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end

            %================
            % Exitwave Update 
            %================
            
            tic
            
%             GPU_phiOLD = gather( sol.GPU.phi );
                     
            [ sol.GPU.phi, ~ ] = RAAR_GPU_arrays_hadamard( sol.GPU.phi,         ...          % RAAR_GPUvec
                                                              sol.GPU.probe,       ...  
                                                              sol.GPU.TFv,         ...
                                                              sol.GPU.ind,         ...
                                                              sol.GPU.sz,          ...
                                                              sol.GPU.Nspos,       ...
                                                              sol.GPU.sqrt_rc,     ...
                                                              sol.GPU.meas_D,      ...
                                                              sol.GPU.meas_Deq0,   ...
                                                              sol.GPU.measLPF,     ...
                                                              sol.GPU.RAAR_beta );

            %========

%             sol.GPU.phi = ER_GPU_arrays_hadamard( sol.GPU.probe, sol.GPU.TFv, sol.GPU.ind, sol.GPU.sz, sol.GPU.Nspos, sol.GPU.sqrt_rc, sol.GPU.meas_D );
                            % ER_GPUvec
                            
 
            %=======
            
%             exwv_diff( gg ) = norm( gather( sol.GPU.phi( : )) - GPU_phiOLD( : ));
     
            %%%%%%%%%%%%%%%
            % Sample Update
            %%%%%%%%%%%%%%%

            if mod( sol.it.exwv, sol.it.sample_update ) == 0

                [ sol.GPU.TFv ] = DMupdate_2Dsample_repmat( sol.GPU.spos_conjP_exwv, ...
                                                            sol.GPU.spos_abs_P_abs2, ...
                                                            sol.GPU.phi,             ...
                                                            sol.GPU.TFv,             ...
                                                            sol.GPU.probe,           ...
                                                            sol.GPU.ind_offset,      ...
                                                            sol.GPU.rc,              ...
                                                            sol.GPU.Nspos );
     

                %====================
                % Sparse Sample Edges
                %====================

                if mod( sol.it.exwv, sol.it.sample_sparsity ) == 0

                    %======

        %             TF  = reshape( sol.GPU.TFv, sol.GPU.samsz );
        %             TF  = lpf_gauss( TF, 0.70 * sol.GPU.samsz );
        %             TF  = sparseFDxy_update_2Dsample( TF, sol.GPU.sparseGPU );      
        %             sol.GPU.TFv = TF( : );   
        %             clear( 'TF' )

                    %======

                    [ sol.GPU.TFv ] = sparseFDxy_update_2Dsample( reshape( sol.GPU.TFv, sol.GPU.samsz ), sol.GPU.sparseGPU );
                    
                    sol.GPU.TFv = sol.GPU.TFv( : );

                end
    
                %======================================================
                % Sample inequality constraints (projection operations)
                %======================================================

                if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

                    sol.GPU.TFv = modulus_limits_project( sol.GPU.TFv, sol.GPU.abs_TF_lim );

        %             sol.GPU.TFv = modulus_limits_project( sol.GPU.TFv, [ 0, 1 ] );
        %             sol.GPU.TFv = modulus_limits_scale( sol.GPU.TFv, sol.GPU.abs_TF_lim );

%                     sol.GPU.TFv = phase_limits_project( sol.GPU.TFv, sol.GPU.phs_TF_lim );
%                     sol.GPU.TFv = phase_limits_scale( sol.GPU.TFv, sol.GPU.phs_TF_lim );

                end
                
                %======================================================
                % Sample mag and phase correlation ( weighted average )
                %======================================================

                if 0 %mod( sol.it.exwv, 50 ) == 0

                    abs_TF = abs( sol.GPU.TFv );
                    abs_TF = abs_TF - min( abs_TF( : ));
                    abs_TF = abs_TF / max( abs_TF( : ));

                    phs_TF = angle( sol.GPU.TFv );
                    phs_TF = phs_TF - min( phs_TF( : ));
                    phs_TF = phs_TF / max( phs_TF( : ));

        %             as = 0.75;
        %             ps = 0.9;
        %             sol.GPU.TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));

                    as = 0.5;
                    sol.GPU.TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.TFv ));

                end  

            end
    
            %%%%%%%%%%%%%%
            % Probe Update 
            %%%%%%%%%%%%%%

            if ( mod( sol.it.exwv, sol.it.probe_update ) == 0 ) && ( sol.it.exwv > sol.it.probe_start )

                %==================================
                % DM/ePIE hybrid style probe update  
                %==================================
                
                TFview = reshape( sol.GPU.TFv( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
                                
                abs2_TFview = abs( TFview ) .^ 2;
                
%                 tmp0 = sum( conj( TFview ) .* sol.GPU.phi - sol.GPU.probe .* abs2_TFview, 4 );
%                 tmp1 = sum( abs2_TFview, 4 );
%                 sol.GPU.probe = sol.GPU.probe + tmp0 ./ ( 1e-7 + tmp1 );

                sol.GPU.probe = sol.GPU.probe + ( sum( conj( TFview ) .* sol.GPU.phi - sol.GPU.probe .* abs2_TFview, 4 ) ) ./ ( 1e-7 + sum( abs2_TFview, 4 ) );
                
                %======================
                % DM style probe update    
                %======================
                
%                 TFview = reshape( sol.GPU.TFv( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
%                 sol.GPU.probe = sum( conj( TFview ) .* sol.GPU.phi, 4 ) ./ sum( 1e-7 + abs( TFview ) .^ 2, 4 );
                
                %==============
                % Probe Support
                %==============
                   
                if ( mod( sol.it.exwv, sol.it.probe_support ) == 0 ) 
                    
                    %========
                    
%                     sol.GPU.probe = sol.GPU.probe .* sol.GPU.probe_support; 
                    
                    %========
                    
                    swparams.blurx     = gpuArray( single( 0.02 )); 
                    swparams.blury     = gpuArray( single( 0.02 ));
                    swparams.sparselvl = gpuArray( single( 0.80 )); 
                    
                     %for pp = 1 : sol.GPU.Nscpm
                     %    
                     %    [ sol.GPU.probe( :, :, pp ), ~, ~ ] = shrinkwrap( sol.GPU.probe( :, :, pp ), swparams ); 
                     %    
                     %end

                    
                    % Support from summed probe mode intensity
                    tmp0          = sum( abs( sol.GPU.probe ) .^ 2, 3 );
                    [ ~, supp ]   = shrinkwrap( tmp0, swparams ); 
                    sol.GPU.probe = sol.GPU.probe .* supp;
                    
%                     [ sol.GPU.probe( :, :, 3 ), supp ] = shrinkwrap( sol.GPU.probe( :, :, 3 ), swparams ); 
%                     sol.GPU.probe = sol.GPU.probe .* supp;
          
                end

                %===================
                % Max abs constraint
                %===================

%                 if 0 %( mod( sol.it.exwv, 1 ) == 0 ) && ~isempty( sol.GPU.scpmmax )
%                     
% %                     tmp2 = reshape( sol.GPU.scpmmax, [ 1, 1, sol.GPU.Nscpm ] );
%                     tmp2 = sol.GPU.scpmmax;
%                     tmp0 = ( abs( sol.GPU.probe ) > tmp2 );
%                     tmp1 = not( tmp0 );
%                         
%                     sol.GPU.probe = sol.GPU.probe .* tmp1 + tmp2 .* exp( 1i * angle( sol.GPU.probe )) .* tmp0;
%                     
%                     clear( 'tmp0', 'tmp1' )
%                     
%                 end
   
                %============================
                % Probe Scaling ( # Photons )
                %============================
                
                if ( mod( sol.it.exwv, sol.it.probe_scaling ) == 0 ) 

                    [ sol.GPU.probe, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.probe, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

                end

                %========================
                % Orthogonalize the SCPMs    % DON'T DO THIS HERE !
                %========================

%                 if ( mod( sol.it.exwv, sol.it.probe_orthog ) == 0 ) 
%                     
%                     [ sol.GPU.probe ] = orthog_modes_eigendecomp( sol.GPU.probe ); 
%                     
%                     
%                     
% %                     tmp0          = abs( sol.GPU.probe );
% %                     [ tmp0 ]      = orthog_modes_eigendecomp( tmp0 ); 
% %                     sol.GPU.probe = tmp0 .* exp( 1i * angle( sol.GPU.probe )); 
%                     
%                 
%                 end


    %             [ scpm ] = compute_scpm_photonocc( sol.GPU.probe );

                   
            end

            %============================
            % Exit wave iteration counter
            %============================

            sol.it.exwv = sol.it.exwv + 1;  
            
            %===========================
            % Fixed Wall time accounting
            %===========================
            
            t = toc;
            sol.total_t = sol.total_t + t;
            
            if sol.total_t > sol.total_t_max
                
                sol.sample.TF = gather( reshape( sol.GPU.TFv, [ sol.sample.sz.sz ]));
                sol.probe.P   = gather( sol.GPU.probe );
                sol.phi       = gather( sol.GPU.phi );
                
                sol = rmfield( sol, 'GPU' );
                
                return
            end
            
        end

        %========================================================
        % Orthogonalize the SCPMs (spatial coherence probe modes)
        %========================================================
                
%         [ sol.GPU.probe ] = orthog_modes_eigendecomp( sol.GPU.probe ); 
        
        %===================
        % GPU to main memory
        %===================

        sol.sample.TF = gather( reshape( sol.GPU.TFv, [ sol.sample.sz.sz ]));
        sol.probe.P   = gather( sol.GPU.probe );
        sol.phi       = gather( sol.GPU.phi );
        
        %========================================================
        % Orthogonalize the SCPMs (spatial coherence probe modes)
        %========================================================
        
        if ( mod( sol.it.exwv - 1, sol.it.probe_orthog ) == 0 )
            
            [ sol.probe.P ] = orthog_modes_eigendecomp( sol.probe.P ); 
            
        end
        
        %======================================================
        % Collect cost function metrics, plot results, clean up
        %======================================================

        if ( mod( sol.it.exwv - 1, sol.it.print_img_results ) == 0 ) 

            [ sol ]       = ptycho2DTPA_collectmetrics( sol, expt, N_datamovements, kk );
            [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );
            
        end

    end
    
    sol = rmfield( sol, 'GPU' );

end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    %==========================================
    
    GPU.RAAR_beta = gpuArray( sol.RAAR.beta );
    
    GPU.batch_N          = gpuArray( sol.spos.batch_N );
    GPU.batch_indxsubset = gpuArray( sol.spos.batch_indxsubset );
    GPU.batch_rs         = gpuArray( sol.spos.batch_rs );
    
    GPU.Nspos = gpuArray( sol.spos.batch_N );
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    %==============
    
    GPU.sz      = gpuArray( sol.sz.sz );
    GPU.rc      = gpuArray( sol.sz.rc );
    GPU.sqrt_rc = gpuArray( sol.sz.sqrt_rc );

    %==============

    GPU.samsz = gpuArray( sol.sample.sz.sz );
    GPU.samrc = gpuArray( sol.sample.sz.rc );

    %==============
    
    GPU.phi = gpuArray( sol.phi );                                                       

    %==============

    % GPU.ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
    
    GPU.ind = gpuArray( sol.spos.frameindx ); 

%     GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = gpuArray( GPU.ind + GPU.ind_offset );
    
    %==============

    % meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = reshape( expt.meas.D( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
    GPU.meas_D = gpuArray( GPU.meas_D );
    
%     GPU.meas_Deq0 = reshape( expt.meas.Deq0( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
%     GPU.meas_Deq0 = gpuArray( GPU.meas_Deq0 );
%     GPU.meas_Deq0 = gpuArray( [] );
    GPU.meas_Deq0 = ( GPU.meas_D == 0 );

    GPU.measLPF = reshape( sol.measLPF, [ GPU.sz, 1, 1 ] );
    GPU.measLPF = gpuArray( sol.measLPF );
    
    %==============

    GPU.probe         = gpuArray( sol.probe.P );
    GPU.fro2TOT       = gpuArray( sol.probe.scpm.fro2TOT );
    GPU.scpmocc       = gpuArray( sol.probe.scpm.occ );
    GPU.probe_support = gpuArray( sol.probe.support );
    GPU.scpmmax       = gpuArray( sol.probe.scpm.max );
    
    %==============
    
    GPU.TFv        = gpuArray( sol.sample.TF( : ));
    GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);

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

%==================================================================================================
