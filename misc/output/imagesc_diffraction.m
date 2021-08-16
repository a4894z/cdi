function imagesc_dffrctn( diffraction_meas, minmax_vals )

if ~exist('minmax_vals','var'),
  minmax_vals(1,1) = 0;
  minmax_vals(1,2) = max(abs(diffraction_meas(:)));
  minmax_vals = log10(1+minmax_vals);
end

imagesc( log10(1+abs(fftshift( diffraction_meas ))), minmax_vals );
daspect([1 1 1])

end

