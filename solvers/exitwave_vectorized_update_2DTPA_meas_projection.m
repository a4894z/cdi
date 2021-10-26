function psi = exitwave_vectorized_update_2DTPA_meas_projection( phi,       ...
                                                                 T0,        ...
                                                                 ind,       ...
                                                                 sz,        ...
                                                                 Nspos,     ...
                                                                 sqrt_rc,   ...
                                                                 meas_D,    ...
                                                                 meas_Deq0, ...
                                                                 measLPF )
    
% exitwave_vectorized_update_2DTPA_meas_projection

    %========================================
    % get sample frames at each scan position
    %========================================
    
    T = reshape( T0( ind ), [ sz, 1, Nspos ]);
 
    %===========================
    % form exitwaves using 2DTPA
    %===========================
    
    psi_PS = T .* phi;

    %==============================================================================
    % measurement constraints assuming Gaussian noise with (ignored) constant stdev
    %==============================================================================

    V = fft( fft( fftshift( fftshift( psi_PS, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
    
    sum_abs2_V = sqrt( sum( abs( V ) .^ 2, 3 ));

    tmp0 = ( V ./ ( 1e-7 + sum_abs2_V ));

    tmp0 = ( meas_D .* tmp0 + V .* meas_Deq0 ) .* measLPF;
%     tmp0 = ( meas_D .* tmp0 + V .* ( meas_D == 0 )) .* measLPF;
    
    psi = fftshift( fftshift( ifft( ifft( tmp0, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    
    %===============================================
    % reflection operator for measurement projection
    %===============================================
    
%     psi = 2 * psi - psi_PS;

end