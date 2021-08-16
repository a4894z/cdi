function [ sol, expt ] = ptycho2DTPA_runiterations( sol, expt, Nit )

%=======================================================
% ???
%=======================================================

sol.spos.batch_indxsubset = 1 : sol.spos.N;

sol.spos.batch_rs  = sol.spos.rs( sol.spos.batch_indxsubset, : );

sol.spos.frameindx = get_indices_2Dframes( sol.spos.batch_rs, sol.sample.sz.sz, sol.sample.vs );

%=======================================================
% form initial exitwaves from current sample/probe modes
%=======================================================

TF = sol.sample.TF( : );
TF = reshape( TF( sol.spos.frameindx ), [ sol.sz.sz, 1, size( sol.spos.frameindx, 2 ) ]);

sol.phi = TF .* sol.probe.P;

%========

[ sol.GPU ] = CPUmem2GPUmem( sol, expt );
% [ sol.GPU ] = CPUmem( sol, expt );

sol.GPU.TF = gpuArray( sol.sample.TF );  

sol.GPU.updateorder = randperm( length( sol.GPU.batch_indxsubset ));

%========

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
                         
                                           
[ sol.sample.TF ] = ePIEupdate_sample( sol.GPU.probe, sol.GPU.TF, sol.GPU.phi, sol.GPU.updateorder, sol.GPU ); 






[ sol, expt ] = update_2DTPAsampleTF( sol, expt );
[ sol, expt ] = update_2DTPAprobe( sol, expt ); 

                                           
%         %=================================
%         %--------- Sample update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.sample_start ) && ( mod( sol.it.exwv, sol.it.sample_update ) == 0 )
%             [ sol, expt ] = update_2DTPAsampleTF( sol, expt ); end
% 
%         %=================================
%         %---------- Probe update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.probe_start ) && ( mod( sol.it.exwv, sol.it.probe_update) == 0 )
%             [ sol, expt ] = update_2DTPAprobe( sol, expt ); end

        
end

%====================================================================================================================================================
    
function [ GPU ] = CPUmem2GPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    %==========================================
    
    GPU.RAAR_beta = gpuArray( sol.RAAR.beta );
    
    GPU.batch_N          = gpuArray( length( sol.spos.batch_indxsubset ) );
    GPU.batch_indxsubset = gpuArray( sol.spos.batch_indxsubset );
    GPU.batch_rs         = gpuArray( sol.spos.batch_rs );
    
    GPU.Nspos = gpuArray( GPU.batch_N );
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

    GPU.ind = uint32( sol.spos.frameindx ); 

%     GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    
    GPU.ind_offset = GPU.ind + GPU.ind_offset;
    
    GPU.ind        = gpuArray( GPU.ind ); 
    GPU.ind_offset = gpuArray( GPU.ind_offset );
    
    %==============

    % meas_D = gpuArray( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = reshape( expt.meas.D( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
    GPU.meas_D = gpuArray( GPU.meas_D );
    
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
    
    GPU.swparams.blurx     = gpuArray( sol.swparams.blurx ); 
    GPU.swparams.blury     = gpuArray( sol.swparams.blury );
    GPU.swparams.sparselvl = gpuArray( sol.swparams.sparselvl ); 
    
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

%====================================================================================================================================================
    
function [ GPU ] = CPUmem( sol, expt )

    %==========================================
    % Move relevant information onto GPU memory
    %==========================================
    
    GPU.RAAR_beta = single( sol.RAAR.beta );
    
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
    
    GPU.phi = single( sol.phi );                                                       

    %==============

    % GPU.ind = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
    
%     GPU.ind = single( sol.spos.frameindx ); 
    GPU.ind = uint32( sol.spos.frameindx ); 
    
%     GPU.ind_offset = double( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
%     GPU.ind_offset = single( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( 0 : GPU.samrc : ( GPU.samrc * ( GPU.Nspos - 1 )));
    GPU.ind_offset = uint32( GPU.ind + GPU.ind_offset );
    
    %==============

    % meas_D = single( permute( repmat( expt.meas.D, 1, 1, 1, GPU.Nscpm ), [ 1, 2, 4, 3 ] ));

    GPU.meas_D = reshape( expt.meas.D( :, :, GPU.batch_indxsubset ), [ GPU.sz, 1, GPU.Nspos ] );
    GPU.meas_D = single( GPU.meas_D );
    
    GPU.meas_Deq0 = ( GPU.meas_D == 0 );

    GPU.measLPF = reshape( sol.measLPF, [ GPU.sz, 1, 1 ] );
    GPU.measLPF = single( sol.measLPF );
    
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
    
    GPU.TFv        = single( sol.sample.TF( : ));
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

    % parallel.gpu.GPUArray.zeros(size(data), 'single');

    GPU.spos_conjP_exwv = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));
    GPU.spos_abs_P_abs2 = single( zeros( [ GPU.samrc, GPU.Nspos ], 'single' ));

    %====================================================
    % Ensure all relevant information is single precision
    %====================================================
    
    % ???

end

%====================================================================================================================================================


            
            
% %=================================================
% %----------- Check if exit wave exists -----------
% %=================================================
%   
%             
%             
%             
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
%------------------------------------------- MAIN LOOP --------------------------------------------
%==================================================================================================
%     
% %     sol.spos.updateorder = randperm( length( sol.spos.indx ), round( 1.0 * length( sol.spos.indx )));
%         
%     for kk = 1 : Nit
% 
%     %     fprintf( '.' );
%         fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
%         if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end
%         
%         %===============================================================
%         %-------- Scan positions we'll update the exit wave over -------
%         %===============================================================
%         
% %         sol.spos.updateorder = randperm( length( sol.spos.indx ), round( 1.0 * length( sol.spos.indx )));
% %         sol.spos.updateorder = sol.spos.indx;
%         
%     %     tmp0 = round( 0.5 * 55 + ( -18 : 18 )); sol.spos.updateorder = tmp0( randperm( length( tmp0 )));
%     
%         if ~isfield( sol.spos, 'updateorder' ) || ( sol.spos.suofai == false ) % same update order for all iterations
% 
%             sol.spos.updateorder = randperm( length( sol.spos.indxsubset ), round( sol.spos.rand_ss * length( sol.spos.indxsubset )));
% %             sol.spos.updateorder = sol.spos.indx;
% 
%         end
% 
%         %=================================
%         %-------- Exit wave update -------
%         %=================================
% 
%         % SPARSE EXIT WAVE OR SPARSE SAMPLE?
%         
%         sol.phi_update = 'RAAR';
% 
%         switch sol.phi_update
% 
%             case 'RAAR'
% 
%                 [ sol.phi ] = RAARupdate_2DTPAexwv_PTFspos_meas( sol.phi, sol, expt );
% 
%             case 'RAAR-ML'
% 
%                 [ sol.phi ] = RAARupdate_2DTPAexwv_PTFspos_PMLmeas( sol.phi, sol, expt );
% 
%             otherwise
% 
%                 [ sol.phi ] = ERupdate_2DTPAexwv_PTFspos_meas( sol, expt );  
%         end
% 
%         %=================================
%         %--------- Sample update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.sample_start ) && ( mod( sol.it.exwv, sol.it.sample_update ) == 0 )
%             [ sol, expt ] = update_2DTPAsampleTF( sol, expt ); end
% 
%         %=================================
%         %---------- Probe update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.probe_start ) && ( mod( sol.it.exwv, sol.it.probe_update) == 0 )
%             [ sol, expt ] = update_2DTPAprobe( sol, expt ); end
% 
%         %===================================================
%         %--------- Scan position correction update ---------
%         %===================================================
% 
% %         if ( sol.it.exwv > sol.it.spos_start  ) && ( mod( sol.it.exwv, sol.it.spos_update ) == 0 )
% % 
% % %             sol.spos.updateorder = randperm( length( sol.spos.indx ), round( 0.2 * length( sol.spos.indx )));     
% % %             sol.spos.update.indx = randperm( length( sol.spos.indx ), round( 0.5 * length( sol.spos.indx )));
% %             sol.spos.update.indx = sol.spos.updateorder;
% %             
% %             for ii = 1 : 1
% % 
% %                 %======================
% % 
% %                 [ sol.spos.rs ] = sposupdateSD_FDexactlinesearch( sol, expt );
% % 
% %                 %======================
% % 
% %                 figure( 666 ); 
% %                 set( gcf, 'Visible', 'off', 'Position', round( [ 1, 1, 0.5 * [ 1920, 1080 ]]))
% %     %                 plot_2Dscan_positions( expt.spos.rs_ALL, [], sol.spos.rs( 1 : 1 : end, : ), sol.spos.indx( 1 : 1 : end ) )
% %                 plot_2Dscan_positions( expt.spos.rs_ALL, [], sol.spos.rs( 1 : 1 : end, : ), [] )
% %                 daspect([1 1 1])
% %                 set( gca, 'xdir', 'reverse' )
% %                 set( gca, 'ydir', 'normal' )
% %                 xlabel('bxh, (+ is left towards door as photons see it)'); ylabel('byv, (+ is up as photons see it)');
% %                 title('positions for scanning the probe on the sample')
% %                 fname = ['scanpos', num2str( [ sol.it.exwv, ii ], '%d-%d' ) ,'.jpg'];
% %                 export_fig( fname, '-r300' ) 
% %                 close all;
% % 
% %             end
% %         end
%           
%         %=================================
%         %-------- Collect metrics --------
%         %=================================
% 
%         if ( mod( sol.it.exwv, sol.it.collect_metrics ) == 0 ) 
%             [ sol ] = ptycho2DTPA_collectmetrics( sol, expt, Nit, kk ); end             % ptycho_collect_metrics()
% 
%         %==================================
%         %----- Print image of results -----
%         %==================================
% 
%         if ( mod( sol.it.exwv, sol.it.print_img_results ) == 0 )
%             [ sol, expt ] = ptycho2DTPA_plotresults( sol, expt ); end               % ptycho_plot_results()
%         
%         %=============================================
%         %----- Exitwave update iteration counter -----
%         %=============================================
% 
%         sol.it.exwv = sol.it.exwv + 1;
% 
%     end
% 
% end

%==================================================================================================
%==================================================================================================
%==================================================================================================
% 
% function [ sol, expt ] = update_2DTPAsampleTF( sol, expt )
% 
%     %=====================================================================
%     %----- e/rPIE or DM style sample update using exitwave and probe -----
%     %=====================================================================
% 
%     %         sol.spos.updateorder = randperm( length( sol.spos.indx ), 1.0 * length( sol.spos.indx ));
% 
%     [ sol.sample.TF ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol ); 
% %     [ sol.sample.TF ] = rPIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol ); 
% %     [ sol.sample.TF ] = mrPIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol, expt );
% 
%     % [ sol.sample.TF ] = DMupdate_sample( sol.probe.P, sol.phi, sol ); 
% 
% 
% %==========================================================
% %-------- Other sample contraints and projections ---------
% %==========================================================
% 
% C1 = sol.it.exwv > sol.it.sample_start;
% 
% %==================
% % Sparse sample constraints 
% 
% C2 = mod( sol.it.exwv, sol.it.sample_sparsity ) == 0;
% 
% if C1 && C2
% 
%     [ sol.sample.TF ] = sparseFDxy_update_2Dsample( sol.sample.TF, sol.sparse ); 
% 
% %             [ sol.sample.TF ] = sparseFDxy_update_2Dsample( sol, sol.sparse_phsTF ); 
% 
% end
%     
% %==================
% % Sample magnitude and phase inequality constraints
% 
% C2 = mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0;
% 
% if C1 && C2
% 
% %         sol.sample.TF = modulus_limits_scale( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     sol.sample.TF = modulus_limits_project( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
% 
% 
% 
%     % !!!!!!!!!!!!!!!!!!!!!!! TRY HISTOGRAM EQUALIZATION CONSTRAINT !!!!!!!!!!!!!!!!!!!!!!!
%     % !!!!!!!!!!!!!!!!!!!!!!! TRY HISTOGRAM EQUALIZATION CONSTRAINT !!!!!!!!!!!!!!!!!!!!!!!
%     % !!!!!!!!!!!!!!!!!!!!!!! TRY HISTOGRAM EQUALIZATION CONSTRAINT !!!!!!!!!!!!!!!!!!!!!!!
%     
% 
% % %         sol.sample.TF = phase_limits_scale( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );
%     sol.sample.TF = phase_limits_project( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );
% 
% %         sol.sample.TF = abs( sol.sample.TF );
% 
% end
% 
% %==================
% % sample mag and phase are correlated?
% 
% % C2 = mod( sol.it.exwv, sol.it.sample_mag_corr ) == 0;
% % 
% % if C1 && C2
% % 
% %     abs_TF = abs( sol.sample.TF );
% %     abs_TF = abs_TF / max( abs_TF( : ));
% % 
% %     phs_TF = angle( sol.sample.TF );
% %     phs_TF = phs_TF - min( phs_TF );
% %     phs_TF = phs_TF / max( phs_TF( : ));
% % 
% %     a = 0.5;
% %     abs_phs_TF = a * abs_TF + ( 1 - a ) * phs_TF;
% % 
% %     sol.sample.TF = abs_TF .* exp( 1i * abs_phs_TF );
% %     sol.sample.TF = abs_phs_TF .* exp( 1i * phs_TF );
% %     sol.sample.TF = abs_phs_TF .* exp( 1i * abs_phs_TF );
% % 
% % end
%     
%     
% end
% 
% % delta_x = ( delta_z / sin( theta ) ) * 0.5 or something
% 
% %==================================================================================================
% %==================================================================================================
% %==================================================================================================
% 
% function [ sol, expt ] = update_2DTPAprobe( sol, expt )
% 
% 
% %         sol.spos.updateorder = randperm( length( sol.spos.indx ), 1.0 * length( sol.spos.indx ));
% 
%     switch sol.probe_update
% 
%         case 'ePIE'
% 
%             [ sol.probe.P ] = ePIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );    
%             
% %             [ sol.probe.P ] = mrPIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );   
%             
%         case 'DM'
%         
%             [ sol.probe.P ] = DMupdate_probemodes( sol.sample.TF, sol.phi, sol );    
%         
%         otherwise
% 
%             warning('Invalid probe update type, defaulting to ePIE.')
%             
%             [ sol.probe.P ] = ePIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );  
%             
%     end
%     
% 
%     
%     
% %==========================================================
% %--------- Other probe contraints and projections ---------
% %==========================================================
% 
% C1 = sol.it.exwv > sol.it.probe_start;
% 
% %==================
% % Probe Support
% 
% C2 = mod( sol.it.exwv, sol.it.probe_support ) == 0;
% 
% if C1 && C2
% 
% %     for pp = 1 : sol.probe.scpm.N
% % 
% % %         swparams.pct_of_max = 0.2;
% % %         swparams.sparselvl = 0.21;
% % %         swparams.blurx = 0.02;
% % %         swparams.blury = 0.02;
% % %         [ sol.probe.P( :, :, pp ), Rthresh ] = shrinkwrap( sol.probe.P( :, :, pp ), swparams );
% % 
% %     
% %     end
%     
%     
%     sol.probe.P = sol.probe.P .* sol.probe.support; 
% 
% end
% 
% %==================
% % Probe Shrinkwrap
% 
% % C2 = mod( sol.it.exwv, sol.it.probe_shrinkwrap ) == 0;
% % 
% % if C1 && C2
% % 
% % %     for pp = 1 : sol.probe.scpm.N
% % % 
% % % %         swparams.pct_of_max = 0.2;
% % % %         swparams.sparselvl = 0.21;
% % % %         swparams.blurx = 0.02;
% % % %         swparams.blury = 0.02;
% % % %         [ sol.probe.P( :, :, pp ), Rthresh ] = shrinkwrap( sol.probe.P( :, :, pp ), swparams );
% % % 
% % %     
% % %     end
% % 
% % end
% 
% %==================
% % Max abs constraint
% 
% % C2 = mod( sol.it.exwv, sol.it.probe_maxvals ) == 0;
% % 
% % if C1 && C2
% %     
% %     for pp = 1 : sol.probe.scpm.N
% % 
% %         tmp0 = ( abs( sol.probe.P( :, :, pp )) >  sol.probe.scpm.max( pp ));
% %         tmp1 = not( tmp0 );
% % 
% %         sol.probe.P( :, :, pp ) = sol.probe.P( :, :, pp ) .* tmp1 + sol.probe.scpm.max( pp ) .* exp( 1i * angle( sol.probe.P( :, :, pp ))) .* tmp0;
% % 
% %     end
% % 
% % end
% 
% %==================
% % Orthogonalize probe modes
% 
% C2 = mod( sol.it.exwv, sol.it.probe_orthog ) == 0;
% 
% if C1 && C2
%     [ sol.probe.P ] = orthog_modes_eigendecomp( sol.probe.P ); end
% 
% %==================
% % LOOK AT WHAT THIS DOES TO MISSING DATA SCALING IN MEASUREMENT PROJECTION
% % Probe scaling, scpm occupancy,
% 
% C2 = mod( sol.it.exwv, sol.it.probe_scaling ) == 0;
% 
% if C1 && C2
%     [ sol.probe.P, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ ); end     
% 
% % 
% % 
% % 
% % 
% 
% 
%     % TRY for tt = 1 : 5, ePIE (T ), ePIE( P ); end
% 
% 
% end

%==================================================================================================
%==================================================================================================
%==================================================================================================








% function [ sol, expt ] = miaoqi202003_mainloop( sol, expt, Nit )
% 
%     %=================================================
%     %----------- Check if exit wave exists -----------
%     %=================================================
% 
%     if ~isfield( sol, 'phi' )
% 
%         sol.phi = zeros( [ sol.sz.sz, sol.probe.scpm.N, sol.spos.N ], 'single' );
% 
%         for ss = 1 : sol.spos.N
% 
%             rs = +sol.spos.rs( ss, : );
% 
%             [ sol.phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +1 * rs, sol.spos.shifttype );
% 
%         end
% 
%     end
% 
%     %=============================================================
%     %------------------------- MAIN LOOP -------------------------
%     %=============================================================
% 
%     for kk = 1 : Nit
% 
%     %     fprintf( '.' );
%         fprintf( [ num2str( sol.it.exwv, '%d ' ), ', ' ] );
%         if mod( sol.it.exwv, 25 ) == 0, fprintf( '\n' ); end
% 
% %         sol.spos.updateorder = randperm( length( sol.spos.indx ), 1.0 * length( sol.spos.indx ));
% 
%     %     tmp0 = round( 0.5 * 55 + ( -18 : 18 )); sol.spos.updateorder = tmp0( randperm( length( tmp0 )));
% 
%         %=================================
%         %-------- Exit wave update -------
%         %=================================
% 
%         % SPARSE EXIT WAVE OR SPARSE SAMPLE?
%         
%         
%         switch sol.phi_update
% 
%             case 'RAAR'
% 
%                 [ sol.phi ] = RAARupdate_2DTPAexwv_PTFspos_meas( sol.phi, sol, expt );
% 
%             case 'RAAR-ML'
% 
%                 [ sol.phi ] = RAARupdate_2DTPAexwv_PTFspos_PMLmeas( sol.phi, sol, expt );
% 
%             otherwise
% 
%                 [ sol.phi ] = ERupdate_2DTPAexwv_PTFspos_meas( sol, expt );  
%         end
% 
%         %=================================
%         %--------- Sample update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.sample_start ) && ( mod( sol.it.exwv, sol.it.sample_updatedate ) == 0 )
%             [ sol, expt ] = update_2DTPAsampleTF( sol, expt ); end
% 
%         %=================================
%         %---------- Probe update ---------
%         %=================================
% 
%         if ( sol.it.exwv > sol.it.probe_start ) && ( mod( sol.it.exwv, sol.it.probe_updatedate) == 0 )
%             [ sol, expt ] = update_2DTPAprobe( sol, expt ); end
% 
%         %===================================================
%         %--------- Scan position correction update ---------
%         %===================================================
% 
%         if ( sol.it.exwv > sol.it.spos_start  ) && ( mod( sol.it.exwv, sol.it.spos_updatedate ) == 0 )
% 
% %             sol.spos.updateorder = randperm( length( sol.spos.indx ), round( 0.2 * length( sol.spos.indx )));
%             
%             sol.spos.update.indx = randperm( length( sol.spos.indx ), round( 0.5 * length( sol.spos.indx )));
%             
%             for ii = 1 : 1
% 
%                 %======================
% 
%                 [ sol.spos.rs ] = sposupdateSD_FDexactlinesearch( sol, expt );
% 
%                 %======================
% 
%                 figure( 666 ); 
%                 set( gcf, 'Visible', 'off', 'Position', round( [ 1, 1, 0.5 * [ 1920, 1080 ]]))
%     %                 plot_2Dscan_positions( expt.spos.rs_ALL, [], sol.spos.rs( 1 : 1 : end, : ), sol.spos.indx( 1 : 1 : end ) )
%                 plot_2Dscan_positions( expt.spos.rs_ALL, [], sol.spos.rs( 1 : 1 : end, : ), [] )
%                 daspect([1 1 1])
%                 set( gca, 'xdir', 'reverse' )
%                 set( gca, 'ydir', 'normal' )
%                 xlabel('bxh, (+ is left towards door as photons see it)'); ylabel('byv, (+ is up as photons see it)');
%                 title('positions for scanning the probe on the sample')
%                 fname = ['scanpos', num2str( [ sol.it.exwv, ii ], '%d-%d' ) ,'.jpg'];
%                 export_fig( fname, '-r300' ) 
%                 close all;
% 
%             end
%         end
%           
%         
%         
%         
%         
%         
%         
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     
%         %=================================
%         %-------- Collect metrics --------
%         %=================================
% 
%         if ( mod( sol.it.exwv, sol.it.collect_metrics ) == 0 ) 
%             [ sol ] = collect_metrics( sol, expt, Nit, kk ); end
% 
%         %==================================
%         %----- Print image of results -----
%         %==================================
% 
%         if ( mod( sol.it.exwv, sol.it.print_img_results ) == 0 )
%             [ sol, expt ] = miaoqi202003_plotresults( sol, expt ); end
% 
%         %=============================================
%         %----- Exitwave update iteration counter -----
%         %=============================================
% 
%         sol.it.exwv = sol.it.exwv + 1;
% 
%     end
% 
% end
% 
% %==================================================================================================
% %==================================================================================================
% %==================================================================================================
% 
% function [ sol, expt ] = update_2DTPAsampleTF( sol, expt )
% 
%     %=====================================================================
%     %----- e/rPIE or DM style sample update using exitwave and probe -----
%     %=====================================================================
% 
%     %         sol.spos.updateorder = randperm( length( sol.spos.indx ), 1.0 * length( sol.spos.indx ));
% 
% %     [ sol.sample.TF ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol ); 
% 
%     % [ sol.sample.TF ] = DMupdate_sample( sol.probe.P, sol.phi, sol ); 
% 
% 
%     [ sol.sample.TF ] = mrPIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol ); 
%     
%     
% 
% %==========================================================
% %-------- Other sample contraints and projections ---------
% %==========================================================
% 
% % 
% % 
% % 
% % 
% 
% C1 = sol.it.exwv > sol.it.sample_start;
% 
% %==================
% % Sparse sample constraints 
% 
% C2 = mod( sol.it.exwv, sol.it.sample_sparsity ) == 0;
% 
% if C1 && C2
% 
%     [ sol.sample.TF ] = sparseFDxy_update_2Dsample( sol.sample, sol.sparse ); 
% 
% %             [ sol.sample.TF ] = sparseFDxy_update_2Dsample( sol, sol.sparse_phsTF ); 
% 
% end
%     
% %==================
% % Sample magnitude and phase inequality constraints
% 
% C2 = mod( sol.it.exwv, sol.it.sample_mag_phs_ineq ) == 0;
% 
% if C1 && C2
% 
% %         sol.sample.TF = modulus_limits_scale( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
%     sol.sample.TF = modulus_limits_project( sol.sample.TF, [ sol.sample.absL, sol.sample.absH ] );
% 
% 
% 
%     % TRY HISTOGRAM EQUALIZATION CONSTRAINT
% 
% 
% 
% %         sol.sample.TF = phase_limits_scale( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );
%     sol.sample.TF = phase_limits_project( sol.sample.TF, [ sol.sample.phsL, sol.sample.phsH ] );
% 
% %         sol.sample.TF = abs( sol.sample.TF );
% 
% end
% 
% % 
% % 
% % 
% % 
% 
%     
%     
% end
% 
% % delta_x = ( delta_z / sin( theta ) ) * 0.5 or something
% 
% %==================================================================================================
% %==================================================================================================
% %==================================================================================================
% 
% function [ sol, expt ] = update_2DTPAprobe( sol, expt )
% 
% 
% %         sol.spos.updateorder = randperm( length( sol.spos.indx ), 1.0 * length( sol.spos.indx ));
% 
%     switch sol.probe_update
% 
%         case 'ePIE'
% 
%             [ sol.probe.P ] = ePIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );       
% 
%         case 'DM'
%         
%             [ sol.probe.P ] = DMupdate_probemodes( sol.sample.TF, sol.phi, sol );    
%         
%         otherwise
% 
%             warning('Invalid probe update type, defaulting to ePIE.')
%             
%             [ sol.probe.P ] = ePIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );  
%             
%     end
%     
% 
%     
%     
% %==========================================================
% %--------- Other probe contraints and projections ---------
% %==========================================================
% 
% % 
% % 
% % 
% % 
% 
% C1 = sol.it.exwv > sol.it.probe_start;
% 
% %==================
% % Probe Support
% 
% C2 = mod( sol.it.exwv, sol.it.probe_orthog ) == 0;
% 
% if C1 && C2
% 
%     for pp = 1 : sol.probe.scpm.N
% 
%         sol.probe.P( :, :, pp ) = sol.probe.P( :, :, pp ) .* sol.probe.support; 
% 
%     end
%     
%     
% %     sol.probe.P = sol.probe.P .* sol.probe.support; 
%     
%     
% 
% end
% 
% %==================    
% % Orthogonalize probe modes
% 
% C2 = mod( sol.it.exwv, sol.it.probe_orthog ) == 0;
% 
% if C1 && C2
%     [ sol.probe.P ] = orthog_modes_eigendecomp( sol.probe.P ); end
% 
% %==================
% % LOOK AT WHAT THIS DOES TO MISSING DATA SCALING IN MEASUREMENT PROJECTION
% % Probe scaling, scpm occupancy,
% 
% C2 = mod( sol.it.exwv, sol.it.probe_scaling ) == 0;
% 
% if C1 && C2
%     [ sol.probe.P, ~, ~ ] = enforce_scpm_fro2TOT_photonocc( sol.probe.P, sol.probe.scpm.fro2TOT, sol.probe.scpm.occ ); end     
% 
% % 
% % 
% % 
% % 
% 
% 
%     % TRY for tt = 1 : 5, ePIE (T ), ePIE( P ); end
% 
% 
% 
% 
% 
%         
% end
% 
% %==================================================================================================
% %==================================================================================================
% %==================================================================================================
% 
% function [ sol ] = collect_metrics( sol, expt, Nit, kk )
% 
%     tmpA = 0;
%     tmpC = 0;
% 
%     for ss = 1 : length( sol.spos.rs )
% 
%         rs = sol.spos.rs( ss, : );
% 
%         [ cview, ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +1 * rs, sol.spos.shifttype );
% 
%         %==============
% 
%         cumsum = 0;
% 
%         for pp = 1 : sol.probe.scpm.N
% 
%             F_phi = fft2( cview( :, :, pp )) / sol.sz.sqrt_rc;
%     %         F_phi = fft2( sol.phi( :, :, pp, ss )) / sol.sz.sqrt_rc;
%             cumsum = cumsum + abs( F_phi ) .^ 2; 
% 
%         end
% 
%         cumsum = cumsum .* expt.meas.SI( ss ).Dneq0;
%         tmpA = tmpA + norm( expt.meas.SI( ss ).D - sqrt( cumsum ), 'fro' );
% 
%         %==============
% 
%         cumsum = 0;
% 
%         for pp = 1 : sol.probe.scpm.N
% 
%             cumsum = cumsum + norm( cview( :, :, pp ) - sol.phi( :, :, pp, ss ), 'fro' );   
% 
%         end
% 
%         tmpC = tmpC + cumsum;
% 
%         %==============
% 
%     end
% 
%     sol.metrics.meas( sol.it.metr )      = tmpA;
%     sol.metrics.exwv_SP( sol.it.metr )   = tmpC;
% 
%     [ scpm ] = compute_scpm_photonocc( sol.probe.P );
%     sol.metrics.scpm_fro2TOT( sol.it.metr ) = scpm.fro2TOT;
%     sol.metrics.scpm_fro2dominant( sol.it.metr ) = scpm.fro2( end );
%     sol.metrics.scpm_fro2others( sol.it.metr ) = scpm.fro2TOT - scpm.fro2( end );
% 
%     sol.metrics.scpmocc_dominant( sol.it.metr ) = sol.metrics.scpm_fro2dominant( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
%     sol.metrics.scpmocc_others( sol.it.metr ) = sol.metrics.scpm_fro2others( sol.it.metr ) / sol.metrics.scpm_fro2TOT( sol.it.metr );
% 
% 
%     % fprintf( [ num2str( [ kk, Nit, sol.it.exwv, sol.metrics.meas( sol.it.metr ), sol.metrics.exwvupdt( sol.it.metr )], ...
%     %             'iteration = %d / %d, iter total = %d, meas metric = %.4f, exit wave update norm = %.4f' ), '\n' ]);
% 
%     fprintf( [ '\n\n', num2str( [ kk, Nit, sol.it.exwv, sol.metrics.meas( sol.it.metr ), sol.metrics.exwv_SP( sol.it.metr )], ...
%                 'iteration = %d / %d, iter total = %d, meas metric = %.2f, exit wave sample probe difference = %.2f' ), '\n\n' ]);
% 
%     sol.it.mtot( sol.it.metr ) = sol.it.exwv;
%     sol.it.metr = sol.it.metr + 1;    
% 
% 
%     figure( 666 ); 
%     set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )     
% 
%     subplot( 2, 1, 1 ); 
%     semilogy( sol.it.mtot, sol.metrics.meas, '-o', 'Linewidth', 2, 'Color', [1.0, 0, 0 ] ); 
%     grid on
%     title('$ \sum_s || \sqrt{D_s} - \left \vert \mathcal{F}[ P ~ T_s ] \right \vert ||_F $','Interpreter','latex');
% 
%     subplot( 2, 1, 2 ); 
%     semilogy( sol.it.mtot, sol.metrics.exwv_SP, '-o', 'Linewidth', 2, 'Color', [0.0, 0, 0.8 ] ); 
%     % title('$ \sum_s || \phi_s - P( \mathbf{r} )  T( \mathbf{r} - \mathbf{r}_s ) ||_F $','Interpreter','latex');
%     grid on
%     title('$ \sum_s || \phi_s -  P ~ T_s ||_F $','Interpreter','latex');
% 
% 
% 
%     % export_fig( num2str( sol.it.exwv, 'meas_metric-%d.jpg' ), '-r90.0' )
%     export_fig( 'metrics_LSmeasPT_LSphiPT.jpg', '-r90.0' )
% 
%     close all;
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
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% % 
% %         sol.meas_error( sol.it.mtot ) = norm( expt.meas.SI( sol.ss ).D - abs( fft2( fftshift( sol.phiB ))) / sol.sz.sqrt_rc  , 'fro' );
% %         sol.it.metr( sol.it.mtot ) = sol.it.exwv;
% %         sol.it.mtot = sol.it.mtot + 1;
% % 
% %         figure( 666 ); 
% %         set( gcf, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
% %         
% % %         subaxis( 2, 3, 1, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 ); 
% %         subplot( 2, 3, [ 1, 4 ] )
% %         plot( sol.it.metr, sol.meas_error ); 
% %         
% % %         subaxis( 2, 3, 2, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
% %         subplot( 2, 3, 2 )
% %         imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( sol.phiB ))); 
% %         xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
% %         title( num2str( sol.it.exwv, '%d' ))
% %         
% % %         subaxis( 2, 3, 3, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.03, 'MB', 0.13 );
% %         subplot( 2, 3, 3 )
% %         imagesc( sol.plot.xaxis, sol.plot.yaxis, log10( 1 + abs( sol.phiT ))); 
% %         xlabel( 'um' ); ylabel( 'um' ); grid on; %daspect([1 1 1])
% %         title( num2str( sol.it.exwv, '%d' ))
% %                 
% %         
% % 
% % % %         subaxis( 2, 3, 4, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
% % %         subplot( 2, 3, 4 )
% % %         imagescHSV( log10( 1 + abs( sol.phiA )) .* exp( 1i * angle( sol.phiA )), sol.plot ); 
% % %         xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
% % 
% % %         subaxis( 2, 3, 5, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
% %         subplot( 2, 3, 5 )
% %         imagescHSV( log10( 1 + abs( sol.phiB )) .* exp( 1i * angle( sol.phiB )), sol.plot ); 
% %         xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
% %         title( num2str( sol.it.exwv, '%d' ))
% %         
% % %         subaxis( 2, 3, 6, 'SpacingVert', 0, 'MR', 0.02, 'ML', 0.04, 'MT', 0.05, 'MB', 0.05 ); 
% %         subplot( 2, 3, 6 )
% %         imagescHSV( log10( 1 + abs( sol.phiT )) .* exp( 1i * angle( sol.phiT )), sol.plot ); 
% %         xlabel( 'um' ); ylabel( 'um' ); grid on; colormap( expt.cm.blj ); %daspect([1 1 1])
% %         
% %         pause( 0.02 )
% %         
% %         export_fig( num2str( [ sol.it.exwv, sol.ss ], 'exitwave_%d_spos-%d.jpg' ), '-r90.0' )
% %         close all;
% % 
% % 
% % 
% % 
% % 
% % 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% end
