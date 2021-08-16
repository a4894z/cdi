function [ sol, expt ] = ptycho2DTPA_runGPUiterations_v2( sol, expt, Nit )


        
    %===========================
    % Compute initial exit waves
    %===========================
    
    sol.spos.frameindx = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
    sol.spos.batch_rs = sol.spos.rs;
    
    sol.spos.frameindx = get_indices_2Dframes( sol.spos.batch_rs, sol.sample.sz.sz, sol.sample.vs );
        
    if ~isfield( sol, 'phi' )

        TF = sol.sample.TF( : );
        TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, sol.spos.N ]);

        sol.phi = TF .* sol.probe.P;

        clear( 'TF' )

    end

    %========================================================
    % Send relevant information from CPU memory to GPU memory
    %========================================================
    
    [ sol.GPU ] = CPUmem2GPUmem( sol, expt );

    for kk = 1 : Nit

        for gg = 1 : 100

            %===============================
            % Feedback for iteration counter
            %===============================
            
            fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
            if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end
            
            %======================================================================================
            %---------------------------------- Exitwave Update -----------------------------------
            %======================================================================================

            sol.GPU.phi = RAAR_GPU_arrays_hadamard( sol.GPU.phi, sol.GPU.probe, sol.GPU.TFv, sol.GPU.ind, sol.GPU.sz, sol.GPU.Nspos, sol.GPU.sqrt_rc, sol.GPU.meas_D );
    %         sol.GPU.phi = ER_GPU_arrays_hadamard( sol.GPU.probe, sol.GPU.TFv, sol.GPU.ind, sol.GPU.sz, sol.GPU.Nspos, sol.GPU.sqrt_rc, sol.GPU.meas_D );

            %======================================================================================
            %---------------------------------- Sample Update -------------------------------------
            %======================================================================================

            if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

                [ sol.GPU.TFv ] = DMupdate_2Dsample_repmat( sol.GPU.spos_conjP_exwv, ...
                                                            sol.GPU.spos_abs_P_abs2, ...
                                                            sol.GPU.phi, ...
                                                            sol.GPU.TFv, ...
                                                            sol.GPU.probe, ...
                                                            sol.GPU.ind_offset, ...
                                                            sol.GPU.rc, ...
                                                            sol.GPU.samsz, ...
                                                            sol.GPU.Nspos );

            end

            %======================================================
            % Sample inequality constraints (projection operations)
            %======================================================

            if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

                sol.GPU.TFv = modulus_limits_project( sol.GPU.TFv, sol.GPU.abs_TF_lim );

    %             sol.GPU.TFv = modulus_limits_project( sol.GPU.TFv, [ 0, 1 ] );
    %             sol.GPU.TFv = modulus_limits_scale( sol.GPU.TFv, sol.GPU.abs_TF_lim );

                sol.GPU.TFv = phase_limits_project( sol.GPU.TFv, sol.GPU.phs_TF_lim );
    %             sol.GPU.TFv = phase_limits_scale( sol.GPU.TFv, sol.GPU.phs_TF_lim );

            end

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
            % Sample mag and phase correlation ( weighted average )
            %======================================================

            if mod( sol.it.exwv, 25e99 ) == 0

                abs_TF = abs( sol.GPU.TFv );
                abs_TF = abs_TF - min( abs_TF( : ));
                abs_TF = abs_TF / max( abs_TF( : ));

                phs_TF = angle( sol.GPU.TFv );
                phs_TF = phs_TF - min( phs_TF( : ));
                phs_TF = phs_TF / max( phs_TF( : ));

    %             as = 0.75;
    %             ps = 0.9;
    %             sol.GPU.TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));

                as = 0.9;
                sol.GPU.TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( sol.GPU.TFv ));

            end

            %======================================================================================
            %------------------------------------ Probe Update ------------------------------------
            %======================================================================================

            if ( mod( sol.it.exwv, sol.it.probe_update ) == 0 ) && ( sol.it.exwv > sol.it.probe_start )

                %===========================
                % DM/ePIE style probe update  
                %===========================
                
                TFview = reshape( sol.GPU.TFv( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
                                
                abs2_TFview = abs( TFview ) .^ 2;
                
                tmp0 = sum( conj( TFview ) .* sol.GPU.phi - sol.GPU.probe .* abs2_TFview, 4 );
                tmp1 = sum( abs2_TFview, 4 );

                sol.GPU.probe = sol.GPU.probe + tmp0 ./ ( 1e-7 + tmp1 );

%                 %======================
%                 % DM style probe update    
%                 %======================
%                 
%                 TFview = reshape( sol.GPU.TFv( sol.GPU.ind ), [ sol.GPU.sz, 1, sol.GPU.Nspos ]);
%                 sol.GPU.probe = sum( conj( TFview ) .* sol.GPU.phi, 4 ) ./ sum( 1e-7 + abs( TFview ) .^ 2, 4 );
                
                %==============
                % Probe Support
                %==============
                   
    %             if ( mod( sol.it.exwv, sol.it.probe_support ) == 0 ) 
    %                 
    %                 sol.GPU.probe = sol.GPU.probe .* sol.GPU.probe_support; 
    %             
    %             end

                %============================
                % Probe Scaling ( # Photons )
                %============================
                
                if ( mod( sol.it.exwv, sol.it.probe_scaling ) == 0 ) 

                    [ sol.GPU.probe, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.GPU.probe, sol.GPU.fro2TOT, sol.GPU.scpmocc ); 

                end

                %=======================
                % Orthogonalize the SCPM
                %=======================

    %             if ( mod( sol.it.exwv, sol.it.probe_orthog ) == 0 ) 
    %                 
    %                 [ sol.GPU.probe ] = orthog_modes_eigendecomp( sol.GPU.probe ); 
    %             
    %             end


    %             [ scpm ] = compute_scpm_photonocc( sol.GPU.probe );


            end

            %============================
            % Exit wave iteration counter
            %============================

            sol.it.exwv = sol.it.exwv + 1;  
            
        end

        %===================
        % GPU to main memory
        %===================

        sol.sample.TF = gather( reshape( sol.GPU.TFv, [ sol.sample.sz.sz ]));
        sol.probe.P   = gather( sol.GPU.probe );
        sol.phi       = gather( sol.GPU.phi );

        [ sol.probe.P ] = orthog_modes_eigendecomp( sol.probe.P ); 

        %======================================================
        % Collect cost function metrics, plot results, clean up
        %======================================================

        [ sol ]       = ptycho2DTPA_collectmetrics( sol, expt, Nit, kk );
        [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );

    end
    
    sol = rmfield( sol, 'GPU' );

end

%==================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )

    GPU.sz      = gpuArray( sol.sz.sz );
    GPU.rc      = gpuArray( sol.sz.rc );
    GPU.sqrt_rc = gpuArray( sol.sz.sqrt_rc );

    %==============

    GPU.samsz = gpuArray( sol.sample.sz.sz );
    GPU.samrc = gpuArray( sol.sample.sz.rc );

    %==============

    GPU.Nspos = gpuArray( sol.spos.N );
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    %==============

    % GPU.ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
    GPU.ind = sol.spos.frameindx;

    % GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    % GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));

    GPU.ind        = gpuArray( GPU.ind );
    GPU.ind_offset = gpuArray( GPU.ind + GPU.ind_offset );

    %==============

    % parallel.gpu.GPUArray.zeros(size(data), 'single');

    GPU.spos_conjP_exwv = gpuArray( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
    GPU.spos_abs_P_abs2 = gpuArray( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));

    %==============

    GPU.probe = gpuArray( sol.probe.P );
    GPU.TFv   = gpuArray( sol.sample.TF( : ));
    GPU.phi   = gpuArray( sol.phi );

    %==============

    GPU.probe_support = gpuArray( sol.probe.support );

    %==============

    GPU.abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
    GPU.phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);

    GPU.fro2TOT = gpuArray( 1.0 * sol.probe.scpm.fro2TOT );
    GPU.scpmocc = gpuArray( sol.probe.scpm.occ );

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

    % meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = gpuArray( expt.meas.D );
    GPU.meas_D = reshape( GPU.meas_D, [ GPU.sz, 1, GPU.Nspos ] );
    
    
    %==============================================================================================
    
    
    
    
    
%     GPU.sz      = single( sol.sz.sz );
%     GPU.rc      = single( sol.sz.rc );
%     GPU.sqrt_rc = single( sol.sz.sqrt_rc );
% 
%     %==============
% 
%     GPU.samsz = single( sol.sample.sz.sz );
%     GPU.samrc = single( sol.sample.sz.rc );
% 
%     %==============
% 
%     GPU.Nspos = single( sol.spos.N );
%     GPU.Nscpm = single( sol.probe.scpm.N );
% 
%     %==============
% 
%     % GPU.ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
%     GPU.ind = sol.spos.frameindx;
% 
%     % GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
% %     GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
% 
%     GPU.ind        = single( GPU.ind );
%     GPU.ind_offset = single( GPU.ind + GPU.ind_offset );
% 
%     %==============
% 
%     % parallel.gpu.GPUArray.zeros(size(data), 'single');
% 
%     GPU.spos_conjP_exwv = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
%     GPU.spos_abs_P_abs2 = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
% 
%     %==============
% 
%     GPU.probe = single( sol.probe.P );
%     GPU.TFv   = single( sol.sample.TF( : ));
%     GPU.phi   = single( sol.phi );
% 
%     %==============
% 
%     GPU.probe_support = single( sol.probe.support );
% 
%     %==============
% 
%     GPU.abs_TF_lim = single( [ sol.sample.absL, sol.sample.absH ]);
%     GPU.phs_TF_lim = single( [ sol.sample.phsL, sol.sample.phsH ]);
% 
%     GPU.fro2TOT = single( 1.0 * sol.probe.scpm.fro2TOT );
%     GPU.scpmocc = single( sol.probe.scpm.occ );
% 
%     %==============
% 
%     GPU.sparseGPU.threshname          = sol.sparse.threshname;
%     GPU.sparseGPU.threshtype          = sol.sparse.threshtype;
%     GPU.sparseGPU.lvl                 = single( sol.sparse.lvl );
%     GPU.sparseGPU.support             = single( sol.sparse.support ); 
%     GPU.sparseGPU.s2DFDxy.fft_scaling = single( sol.sparse.s2DFDxy.fft_scaling );
%     GPU.sparseGPU.s2DFDxy.qLPF        = single( sol.sparse.s2DFDxy.qLPF );
%     GPU.sparseGPU.s2DFDxy.Dx_fft      = single( sol.sparse.s2DFDxy.Dx_fft );
%     GPU.sparseGPU.s2DFDxy.Dy_fft      = single( sol.sparse.s2DFDxy.Dy_fft );
%     GPU.sparseGPU.s2DFDxy.Dx_fft_conj = single( sol.sparse.s2DFDxy.Dx_fft_conj );
%     GPU.sparseGPU.s2DFDxy.Dy_fft_conj = single( sol.sparse.s2DFDxy.Dy_fft_conj );
%     GPU.sparseGPU.s2DFDxy.sqrt_rc     = single( sol.sparse.s2DFDxy.sqrt_rc );
% 
%     %==============
% 
%     % meas_D = single( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));
% 
%     GPU.meas_D = single( expt.meas.D );
%     GPU.meas_D = reshape( GPU.meas_D, [ sol.sz.sz, 1, sol.spos.N ] );
    
    
    
    
    
end

%==================================================================================================










































%{

%==================================================================================================

sz      = gpuArray( sol.sz.sz );
rc      = gpuArray( sol.sz.rc );
sqrt_rc = gpuArray( sol.sz.sqrt_rc );

%==============

samsz = gpuArray( sol.sample.sz.sz );
samrc = gpuArray( sol.sample.sz.rc );

%==============

Nspos = gpuArray( sol.spos.N );
Nscpm = gpuArray( sol.probe.scpm.N );

%==============

% ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
ind = sol.spos.frameindx;

% ind_offset = double( 0 : samrc : ( samrc * ( Nspos - 1 )));
% ind_offset = single( 0 : samrc : ( samrc * ( Nspos - 1 )));
ind_offset = uint32( 0 : samrc : ( samrc * ( Nspos - 1 )));

ind        = gpuArray( ind );
ind_offset = gpuArray( ind + ind_offset );

%==============

% parallel.gpu.GPUArray.zeros(size(data), 'single');

spos_conjP_exwv = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
spos_abs_P_abs2 = gpuArray( zeros( [ samrc, Nspos ], 'single' ));

%==============

probe = gpuArray( sol.probe.P );
TFv   = gpuArray( sol.sample.TF( : ));
phi   = gpuArray( sol.phi );

%==============

probe_support = gpuArray( sol.probe.support );

%==============

abs_TF_lim = gpuArray( [ sol.sample.absL, sol.sample.absH ]);
phs_TF_lim = gpuArray( [ sol.sample.phsL, sol.sample.phsH ]);

fro2TOT = gpuArray( 1.0 * sol.probe.scpm.fro2TOT );
scpmocc = gpuArray( sol.probe.scpm.occ );

%==============

sparseGPU.threshname          = sol.sparse.threshname;
sparseGPU.threshtype          = sol.sparse.threshtype;
sparseGPU.lvl                 = gpuArray( sol.sparse.lvl );
sparseGPU.support             = gpuArray( sol.sparse.support ); 
sparseGPU.s2DFDxy.fft_scaling = gpuArray( sol.sparse.s2DFDxy.fft_scaling );
sparseGPU.s2DFDxy.qLPF        = gpuArray( sol.sparse.s2DFDxy.qLPF );
sparseGPU.s2DFDxy.Dx_fft      = gpuArray( sol.sparse.s2DFDxy.Dx_fft );
sparseGPU.s2DFDxy.Dy_fft      = gpuArray( sol.sparse.s2DFDxy.Dy_fft );
sparseGPU.s2DFDxy.Dx_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dx_fft_conj );
sparseGPU.s2DFDxy.Dy_fft_conj = gpuArray( sol.sparse.s2DFDxy.Dy_fft_conj );
sparseGPU.s2DFDxy.sqrt_rc     = gpuArray( sol.sparse.s2DFDxy.sqrt_rc );

%==============

% meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, Nscpm ), [ 1, 2, 4, 3 ] ));

meas_D = gpuArray( expt.meas.D );
meas_D = reshape( meas_D, [ sol.sz.sz, 1, sol.spos.N ] );

%==============================================================================================
%------------------------------------------- MAIN LOOP ----------------------------------------
%==============================================================================================

for kk = 1 : Nit

    for gg = 1 : 100

        %===============================
        % Feedback for iteration counter
        %===============================

        fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
        if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end
        
        %=================
        % Exitwave updates
        %=================

        phi = RAAR_GPU_arrays_hadamard( phi, probe, TFv, ind, sz, Nspos, sqrt_rc, meas_D );
%         phi = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, sqrt_rc, meas_D );
        
%         phi = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, Nscpm, meas_D );

        %================================================
        % Sample update using current probe and exit wave
        %================================================
        
        if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

            [ TFv ] = DMupdate_2Dsample_repmat( spos_conjP_exwv, ...
                                                spos_abs_P_abs2, ...
                                                phi, ...
                                                TFv, ...
                                                probe, ...
                                                ind_offset, ...
                                                rc, ...
                                                samsz, ...
                                                Nspos );

        end
        
        %======================================================
        % Sample inequality constraints (projection operations)
        %======================================================
        
        if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

            TFv = modulus_limits_project( TFv, abs_TF_lim );

%             TFv = modulus_limits_project( TFv, [ 0, 1 ] );
%             TFv = modulus_limits_scale( TFv, abs_TF_lim );

            TFv = phase_limits_project( TFv, phs_TF_lim );
%             TFv = phase_limits_scale( TFv, phs_TF_lim );

        end

        %====================
        % Sparse Sample Edges
        %====================
        
        if mod( sol.it.exwv, sol.it.sample_sparsity ) == 0
        
            %======

%             TF  = reshape( TFv, samsz );
%             TF  = lpf_gauss( TF, 0.70 * samsz );
%             TF  = sparseFDxy_update_2Dsample( TF, sparseGPU );      
%             TFv = TF( : );   
%             clear( 'TF' )
            
            %======

            [ TFv ] = sparseFDxy_update_2Dsample( reshape( TFv, samsz ), sparseGPU );
            TFv = TFv( : );
  
        end

        %======================================================
        % sample mag and phase correlation ( weighted average )
        %======================================================
        
        if mod( sol.it.exwv, 25e99 ) == 0

            abs_TF = abs( TFv );
            abs_TF = abs_TF - min( abs_TF( : ));
            abs_TF = abs_TF / max( abs_TF( : ));

            phs_TF = angle( TFv );
            phs_TF = phs_TF - min( phs_TF( : ));
            phs_TF = phs_TF / max( phs_TF( : ));

%             as = 0.75;
%             ps = 0.9;
%             TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * 2 * pi * (  ps * abs_TF + ( 1 - ps ) * phs_TF ));
    
            as = 0.9;
            TFv = ( as * abs_TF + ( 1 - as ) * phs_TF ) .* exp( 1i * angle( TFv ));

        end
        
        %=============
        % Probe Update
        %=============

        if ( mod( sol.it.exwv, sol.it.probe_update ) == 0 ) && ( sol.it.exwv > sol.it.probe_start )
            
%             phi = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, sqrt_rc, meas_D );
            
            % !!!!!!!!!!!!!!!!!!!!! MAKE THIS SIMILAR TO SAMPLE DM/ePIE update
            
            TFview = reshape( TFv( ind ), [ sz, Nspos ]);
            tmp0   = permute( repmat( conj( TFview ), [ 1, 1, 1, Nscpm ] ), [ 1, 2, 4, 3 ] );            % ###################################
            probe  = sum( tmp0 .* phi, 4 ) ./ sum( 1e-7 + abs( tmp0 ) .^ 2, 4 );
            clear( 'tmp0', 'TFview' )

% 
%             if ( mod( sol.it.exwv, sol.it.probe_support ) == 0 ) 
%                 
%                 probe = probe .* probe_support; 
%             
%             end

            if ( mod( sol.it.exwv, sol.it.probe_scaling ) == 0 ) 
                
                [ probe, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( probe, fro2TOT, scpmocc ); 
            
            end
            
%             if ( mod( sol.it.exwv, sol.it.probe_orthog ) == 0 ) 
%                 
%                 [ probe ] = orthog_modes_eigendecomp( probe ); 
%             
%             end
            
            
%             [ scpm ] = compute_scpm_photonocc( probe );
            
            
        end

        %============================
        % Exit wave iteration counter
        %============================
        
        sol.it.exwv = sol.it.exwv + 1;  

    end

    %===================
    % GPU to main memory
    %===================
    
    sol.sample.TF = gather( reshape( TFv, [ sol.sample.sz.sz ]));
    sol.probe.P   = gather( probe );
    sol.phi       = gather( phi );
    
    [ sol.probe.P ] = orthog_modes_eigendecomp( sol.probe.P ); 
    
    %============================================
    % Collect cost function metrics; plot results
    %============================================
    
    [ sol ]       = ptycho2DTPA_collectmetrics( sol, expt, Nit, kk );
    [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );

end

%==================================================================================================

%}