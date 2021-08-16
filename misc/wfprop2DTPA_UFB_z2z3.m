function [ phiz3 ] = wfprop2DTPA_UFB_z2z3( phiz2, expt )

% wfpropUFB_z2z3        % pass options for z2z3fresnel, fft2
% wfpropFB_z2z3         % pass options for z2z3fresnel, z2z1z3fresnel, fft2

%==================================================================================================

% propagate all the exitwaves at sample plane z2 ( generated 
% from the SCPMs ) to the measurement plane z3:

%==================================================================================================

% GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP
% GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP
% GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP GET RID OF FOR LOOP

phiz3 = zeros( expt.sz.r, expt.sz.c, expt.probe.scpm.N, 'single' );

for pp = 1 : expt.probe.scpm.N
    
    % using simple fft ( Fraunhofer approx ):
    phiz3( :, :, pp ) = fft2( fftshift( phiz2( :, :, pp ))) / expt.sz.sqrt_rc;

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

