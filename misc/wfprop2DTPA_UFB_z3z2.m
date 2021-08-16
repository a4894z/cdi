function [ phiz2 ] = wfprop2DTPA_UFB_z3z2( phiz3, expt )


% propagate all the exitwaves at sample plane z2 ( generated 
% from the SCPMs ) to the measurement plane z3:

%==================================================================================================

phiz2 = zeros( expt.sz.r, expt.sz.c, expt.probe.scpm.N, 'single' );

%===============================

for pp = 1 : expt.probe.scpm.N
    
    % using simple fft ( Fraunhofer approx ):
    phiz2( :, :, pp ) = fftshift( ifft2( phiz3( :, :, pp ))) * expt.sz.sqrt_rc;

end

%===============================

% wfp.lambda = expt.lambda; 
% wfp.zif = expt.csys.z23; 
% wfp.zi = expt.csys.z2;
% wfp.zf = expt.csys.z3;
% wfp.intcurv = true;
% wfp.extcurv = false;
% wfp.dir = 'forward';
%
% for pp = 1 : expt.probe.scpm.N
% 
%     phiz3( :, :, pp ) = fresnelpropdirect( phiz2, wfp );
% 
% end

%==================================================================================================

