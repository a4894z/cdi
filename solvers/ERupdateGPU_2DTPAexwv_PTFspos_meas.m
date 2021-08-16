function [ sample, probe, phi ] = ERupdateGPU_2DTPAexwv_PTFspos_meas( sol, expt )

%==================================================================================================

    [ ind ] = get_indices_2Dframes( sol.spos.rs, sol.sample.sz.sz, sol.sample.vs );
        
    probe = gpuArray( sol.probe.P );
    TFv   = gpuArray( sol.sample.TF( : ));
    ind   = gpuArray( ind );
    
    sz    = gpuArray( sol.sz.sz );
%     rc    = sz( 1 ) * sz( 2 );

    Nspos = gpuArray( sol.spos.N );
    Nscpm = gpuArray( sol.probe.scpm.N );
    
    meas_D = permute( repmat( gpuArray( expt.meas.D ), 1, 1, 1, Nscpm ), [ 1, 2, 4, 3 ] );
    
    
    
    
    samsz = gpuArray( sol.sample.sz.sz );
%     samrc = samsz( 1 ) * samsz( 2 );
    samrc = gpuArray( sol.sample.sz.rc );
    


    spos_conjP_exwv = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
    spos_abs_P_abs2 = gpuArray( zeros( [ samrc, Nspos ], 'single' ));
    
    ind_offset = uint32( 0 : samrc : ( samrc * ( Nspos - 1 )));
%     ind_offset = 0 : samrc : ( samrc * ( spos_N - 1 ));

%==================================================================================================


    % exit wave update
    phi = ER_GPU_arrays_hadamard( probe, TFv, ind, sz, Nspos, Nscpm, meas_D );


    % sample update, order doesn't matter
    [ sample ] = DMupdate_2Dsample_repmat( spos_conjP_exwv, ...
                                           spos_abs_P_abs2, ...
                                           phiVg, ...
                                           probe, ...
                                           ind + ind_offset, ...
                                           sz, ...
                                           samsz, ...
                                           Nspos );


    % probe update, order doesn't matter
    tmp0   = sample( : );
    TFview = reshape( tmp0( ind ), [ sz, Nspos ]);
    tmp0   = permute( repmat( conj( TFview ), [ 1, 1, 1, Nscpm ] ), [ 1, 2, 4, 3 ] );
    probe  = sum( tmp0 .* phi, 4 ) ./ sum( 1e-7 + abs( tmp0 ) .^ 2, 4 );



    clear( 'probe', 'TFv', 'ind', 'sz', 'Nspos', 'Nscpm', 'meas_D' )
    
%==================================================================================================