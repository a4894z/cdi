function [ psi ] = RAAR_GPU_arrays_hadamard( psi, phi, T, ind, sz, Nspos, sqrt_rc, meas_D, beamstop, measLPF, RAAR_beta )
     
    %========
    
    T = reshape( T( ind ), [ sz, Nspos ]);
    T = reshape( T, [ sz, 1, Nspos ] );

    tmp0 = T .* phi;
    psi_RS = 2 * tmp0 - psi;
     
    %========

    V = fft( fft( fftshift( fftshift( psi_RS, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
    
    tmp0 = sqrt( sum( abs( V ) .^ 2, 3 ));

    tmp0 = ( V ./ ( 1e-7 + tmp0 ));

    tmp0 = meas_D .* tmp0 + V .* beamstop;

    tmp0 = tmp0 .* measLPF;

    tmp0 = fftshift( fftshift( ifft( ifft( tmp0, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    
    psi_RM_RS = 2 * tmp0 - psi_RS;

    %========
    
    V = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
    
    tmp0 = sqrt( sum( abs( V ) .^ 2, 3 ));

    tmp0 = ( V ./ ( 1e-7 + tmp0 ));

    tmp0 = meas_D .* tmp0 + V .* beamstop;
    
    tmp0 = tmp0 .* measLPF;
    
    psi_PM = fftshift( fftshift( ifft( ifft( tmp0, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    
    %========
    
    psi = RAAR_beta * 0.5 * ( psi_RM_RS + psi ) + ( 1 - RAAR_beta ) * psi_PM;
    
end


