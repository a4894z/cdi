function [ psi, metrics ] = exitwave_vectorized_update_2DTPA_meas_projection( phi,          ...
                                                                              T0,           ...
                                                                              ind,          ...
                                                                              sz,           ...
                                                                              Nspos,        ...
                                                                              sqrt_rc,      ...
                                                                              sqrt_I_m,     ...
                                                                              sqrt_I_m_eq0, ...
                                                                              measLPF,      ...
                                                                              collect_metrics )
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
metrics = [];

%====================================================================
% get sample frames at each scan position, form exitwaves using 2DTPA
%====================================================================

psi = reshape( T0( ind ), [ sz, 1, Nspos ]) .* phi;

%==============================================================================
% measurement constraints assuming Gaussian noise with (ignored) constant stdev
%==============================================================================

psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;

%========

if collect_metrics 
    
    I_e      = sum( abs( psi ) .^ 2, 3 );
    sqrt_I_e = sqrt( I_e );
    I_m      = sqrt_I_m .^ 2;
    meas_mask = not( sqrt_I_m_eq0 );
    
    metrics.gauss_intensity = sum( sum( sum( meas_mask .* abs( I_m - I_e ) .^ 2           )));
    metrics.gauss_magnitude = sum( sum( sum( meas_mask .* abs( sqrt_I_m - sqrt_I_e ) .^ 2 )));
    metrics.poiss           = sum( sum( sum( meas_mask .* ( I_e - I_m .* log( I_e ))      )));
    
else
    
    sqrt_I_e = sqrt( sum( abs( psi ) .^ 2, 3 ));

end

%========

psi = ( sqrt_I_m .* ( psi ./ ( 1e-7 + sqrt_I_e )) + psi .* sqrt_I_m_eq0 ) .* measLPF;

psi = fftshift( fftshift( ifft( ifft( psi, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    










                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
% metrics = [];
% 
% %====================================================================
% % get sample frames at each scan position, form exitwaves using 2DTPA
% %====================================================================
% 
% psi = reshape( T0( ind ), [ sz, 1, Nspos ]) .* phi;
% 
% %==============================================================================
% % measurement constraints assuming Gaussian noise with (ignored) constant stdev
% %==============================================================================
% 
% psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
% 
% sum_abs2_Psi = sqrt( sum( abs( psi ) .^ 2, 3 ));
% 
% %========
% 
% if collect_metrics 
%     
%     metrics = sum( sum( sum( abs( meas_D - sum_abs2_Psi .* not( meas_Deq0 ) ) .^ 2 )));
%     
% %     metrics.gauss = sum( sum( sum( abs( sqrt( I_m ) - sqrt( I_e ) .* not( meas_Deq0 ) ) .^ 2 )));
% %     metrics.poiss = sum( sum( sum( I_e -  ( meas_D .^ 2 ) .* log( I_e ))));
% 
% end
% 
% %========
% 
% psi = ( meas_D .* ( psi ./ ( 1e-7 + sum_abs2_Psi )) + psi .* meas_Deq0 ) .* measLPF;
% 
% psi = fftshift( fftshift( ifft( ifft( psi, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;
    
end




