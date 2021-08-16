function [ sol, expt ] = ptycho2DTPA_runGPUiterations( sol, expt, Nit )

% MOD THIS SO WE CAN USE A SUBSET OF SPOS FOR UPDATE...LOAD ALL MEAS ONTO GPU, BUT ONLY USE REPMAT ETC ON SUBSET
 
% if ~isfield( sol, 'phi' )
% 
%     sol.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N, sol.spos.N ], 'single' );
% 
%     for ss = 1 : sol.spos.N
% 
%         rs = +sol.spos.rs( ss, : );
%                                             % pass a subset of spos and loop inside here:
%         [ sol.phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +1 * rs, sol.spos.shifttype );
% 
%     end
% 
% end

%==================================================================================================

sz    = gpuArray( sol.sz.sz );
rc    = gpuArray( sol.sz.rc );

%==============

samsz = gpuArray( sol.sample.sz.sz );
samrc = gpuArray( sol.sample.sz.rc );

%==============

Nspos = gpuArray( sol.spos.N );
Nscpm = gpuArray( sol.probe.scpm.N );

%==============

ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );

% ind_offset = double( 0 : samrc : ( samrc * ( Nspos - 1 )));
% ind_offset = single( 0 : samrc : ( samrc * ( Nspos - 1 )));
ind_offset = uint32( 0 : samrc : ( samrc * ( Nspos - 1 )));

ind        = gpuArray( ind );
ind_offset = gpuArray( ind_offset );

%==============

% parallel.gpu.GPUArray.zeros(size(data), 'single');

spos_conjP_exwv = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
spos_abs_P_abs2 = gpuArray( zeros( [ samrc, Nspos ], 'single' ));

%==============

probe = gpuArray( sol.probe.P );
TFv   = gpuArray( sol.sample.TF( : ));

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

meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, Nscpm ), [ 1, 2, 4, 3 ] ));

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
        
        phiVg = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, Nscpm, meas_D );

        %================================================
        % Sample update using current probe and exit wave
        %================================================
        
        if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

            [ TFv ] = DMupdate_2Dsample_repmat( spos_conjP_exwv, ...
                                                spos_abs_P_abs2, ...
                                                phiVg, ...
                                                TFv, ...
                                                probe, ...
                                                ind + ind_offset, ...
                                                rc, ...
                                                samsz, ...
                                                Nspos );

        end

        %======================================================
        % Sample inequality constraints (projection operations)
        %======================================================
        
        if mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0

            TFv = modulus_limits_project( TFv, abs_TF_lim );
        %     TFv = modulus_limits_scale( TFv, abs_TF_lim );

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

            % tmp0   = sample4b( : );
            % tmp0   = sample4b;
            TFview = reshape( TFv( ind ), [ sz, Nspos ]);
            tmp0   = permute( repmat( conj( TFview ), [ 1, 1, 1, Nscpm ] ), [ 1, 2, 4, 3 ] );
            probe  = sum( tmp0 .* phiVg, 4 ) ./ sum( 1e-7 + abs( tmp0 ) .^ 2, 4 );
            clear( 'tmp0', 'TFview' )

            if ( mod( sol.it.exwv, sol.it.probe_scaling ) == 0 ) 
                [ probe, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( probe, fro2TOT, scpmocc ); end
            
            if ( mod( sol.it.exwv, sol.it.probe_orthog ) == 0 ) 
                [ probe ] = orthog_modes_eigendecomp( probe ); end

            if ( mod( sol.it.exwv, sol.it.probe_support ) == 0 ) 
                probe = probe .* probe_support; end
            
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

    %============================================
    % Collect cost function metrics; plot results
    %============================================
    
    [ sol ]       = ptycho2DTPA_collectmetrics( sol, expt, Nit, kk );
    [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt );

end
    
    
    