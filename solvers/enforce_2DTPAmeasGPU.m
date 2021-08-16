function V = enforce_2DTPAmeasGPU( phi, meas, measLPF, sol )




V = fft( fftshift( fft( fftshift( phi, 1 ), [], 1 ), 2 ), [], 2 ) / sol.sz.sqrt_rc;
tmp0 = sqrt( sum( abs( V ) .^ 2, 3 ));