function [ meas ] = make_measurements_2DTPA_FB( expt )

% define the exitwave at sample plane z2 (assuming projection 
% approximation) and propagate from z2 to measurement plane z3 
% using the "direct" method for long range propagation

%================================================
%------------ sample out measurement ------------
%================================================
 
meas.SO = fresnel_prop_direct( expt.probe.P, expt.z23, expt.L2, expt.L3, expt.lambda, expt.sz );
meas.SO = abs( fftshift( meas.SO ));

% tmp2 = abs( fft2( expt.probe.P )) / expt.sz.sqrt_rc;

%================================================
%------------ sample in measurements ------------
%================================================

for ii = 1 : expt.spos.N

    [ tmp1 ] = projapprox_2Dexitwave( expt.probe.P, expt.sample.T, expt.spos.rj( ii, : ));
    tmp1 = fresnel_prop_direct( tmp1, expt.z23, expt.L2, expt.L3, expt.lambda, expt.sz );

    meas.SI( :, :, ii ) = abs( fftshift( tmp1 ));
        
%     tmp1 = abs( fft2( fftshift( exitwave ))) / expt.sz.sqrt_rc;
    
%     figure;
%     imagescHSV( fftshift( log10(1+abs(meas.SI( :, :, ii )))  .* exp(1i*angle(meas.SI( :, :, ii )))))
%     figure;
%     imagescHSV( fftshift( log10(1+abs(tmp1))  .* exp(1i*angle(tmp1))))
    
end

clear( 'tmp1', 'ii' )

%==================================================================================================
%------------------------- introduce noise to the measurements ------------------------------------
%==================================================================================================

%{
% measurments just defined are noise free, now introduce poisson noise + detector background noise, 
% i.e. a noise model like: In = Poiss[ I ] + Unif[ ]

meas.noisy = true;

% total number of exposures for sample out measurement:
meas.NexpSO = 10;

% total number of exposures for sample in measurement:
meas.NexpSI = 10;

% magnitude --> intensities
measI_SO = single( meas.SO .^ 2 );
measI_SI = single( meas.SI .^ 2 );

% define max background value allowed in unif random distribution:
max_BGSO = 5e-5 * max( measI_SO( : ));
max_BGSI = max_BGSO;                        % should be the same as we generally use same exposure time for sample in/out
% max_BGSI = 5e-5 * max( measI_SI( : ));

%================================================
%-------- for the sample OUT measurement --------
%================================================

tmpSO = zeros( expt.sz.r, expt.sz.c, 'single' );

for ii = 1 : meas.NexpSO

    % TODO: add in an effect where the wavefield before the focusing optic has some variation shot to shot
    % e.g. a blurred random array with variations of +/- 5% or so?
    % or a gaussian of some fwhm that can "jitter" spatially wrt optic position

    % poisson noise in intensity measurements:
    tmp2 = random( 'poisson', measI_SO );

    % constant noisy background ( uniform random distributed ):
    tmp1 = random( 'unif', 0, max_BGSO, expt.sz.r, expt.sz.c );

    % accumulate single exposure intensity measurements:
    tmpSO = tmpSO + tmp2 + tmp1;

end

meas.SO = sqrt( tmpSO / meas.NexpSO );

%================================================
%-------- for the sample IN measurements --------
%================================================

tmpSI = zeros( expt.sz.r, expt.sz.c, 'single' );

for jj = 1 : expt.spos.N
    
    for ii = 1 : meas.NexpSI  

        % TODO: add in an effect where the wavefield before the focusing optic has some variation shot to shot
        % e.g. a blurred random array with variations of +/- 5% or so?
        % or a gaussian of some fwhm that can "jitter" spatially wrt optic position

        % poisson noise in intensity measurements:
        tmp1 = random( 'poisson', measI_SI( :, :, jj ) );

        % constant noisy background ( uniform random distributed ):
        tmp2 = random( 'unif', 0, max_BGSI, expt.sz.r, expt.sz.c );
 
        % accumulate single exposure intensity measurements:
        tmpSI = tmpSI + tmp1 + tmp2;

    end
    
    meas.SI( :, :, jj ) = sqrt( tmpSI / meas.NexpSI );
    tmpSI = zeros( expt.sz.r, expt.sz.c, 'single' );
    
end



% constant background subtraction on intensity:

bg_sub = 70.0;
meas.SI = abs( meas.SI ) .^ 2 - bg_sub;
meas.SI( meas.SI < 0 ) = 0;
meas.SI = sqrt( meas.SI + 0 * (meas.SI ~= 0) * bg_sub );

bg_sub = 70.0;
meas.SO = abs( meas.SO ) .^ 2 - bg_sub;
meas.SO( meas.SO < 0 ) = 0;
meas.SO = sqrt( meas.SO + 0 * (meas.SO ~= 0) * bg_sub );




%================================================
%------------------ clean up --------------------
%================================================


clear( 'tmp*', 'measI_S*', 'max_BGS*', 'ii', 'jj', 'bg_sub' )

%}