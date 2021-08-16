function phiN = ER_CPU_arrays_hadamard( sol, expt )


    [ ind ] = get_indices_2Dframes( sol.spos.rs, expt.sample.sz.sz, expt.sample.vs );

    TFview = sol.sample.TF( : );
    TFview = TFview( ind );
    TFview = reshape( TFview, [ expt.sz.sz, expt.spos.N ]);

    probe  = repmat( sol.probe.P, 1, 1, 1, sol.spos.N  );
    TFview = permute( repmat( TFview, [ 1, 1, 1, sol.probe.scpm.N ] ), [ 1, 2, 4, 3 ] );

    phiN = TFview .* probe;


%     clear( 'TFview', 'probe' )

    %======================================================

    % meas_Deq0 = ( meas_D == 0 );

    meas_D = repmat( expt.meas.D, 1, 1, 1, sol.probe.scpm.N );
    meas_D = permute( meas_D, [ 1, 2, 4, 3 ] );

    V = fft( fft( fftshift( fftshift( phiN, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;
    % V = fft( fftshift( fft( fftshift( phiN, 1 ), [], 1 ), 2 ), [], 2 ) / sol.sz.sqrt_rc;

    sum_abs2_V = repmat( sqrt( sum( abs( V ) .^ 2, 3 )), [ 1, 1, sol.probe.scpm.N, 1 ] );

    tmp1 = meas_D .* ( V ./ ( 1e-7 + sum_abs2_V )) + V .* ( meas_D == 0 );

    phiN = fftshift( fftshift( ifft( ifft( tmp1, [], 1 ), [], 2 ), 1 ), 2 ) * sol.sz.sqrt_rc;

%     clear( 'tmp1', 'meas_D', 'V' , 'sum_abs2_V' )

end