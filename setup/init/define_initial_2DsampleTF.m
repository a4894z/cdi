function [ sol, expt ] = define_initial_2DsampleTF( sol, expt )

%=====================================
% load previous sample here if desired
%=====================================

% [ sol, expt ] = load_previous_2DsampleTF( sol, expt );
% return

%=========================
% zero pad existing sample
%=========================

% sol.sample.T = padarray( sol.sample.T, [ 400, 0 ], 1 );
% 
% sol.sample.sz.r = size( sol.sample.T, 1 );
% sol.sample.sz.c = size( sol.sample.T, 2 );
% 
% sol.sample.sz.sz      = [ sol.sample.sz.r, sol.sample.sz.c ];
% sol.sample.sz.rc      = sol.sample.sz.r * sol.sample.sz.c;
% sol.sample.sz.sqrt_rc = sqrt( sol.sample.sz.rc );
% 
% sol.sample.vs.r = single( round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r ))));
% sol.sample.vs.c = single( round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c ))));
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reinit sample to blurred random numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================
% array size for the sample
%==========================

%{

sol.sample.sz.r = round( 1.00 * round( max( expt.spos.rs( :, 1 )) - min( expt.spos.rs( :, 1 )) + sol.sz.r ) );
sol.sample.sz.c = round( 1.00 * round( max( expt.spos.rs( :, 2 )) - min( expt.spos.rs( :, 2 )) + sol.sz.c ) );
% sol.sample.sz.r = round( 1.00 * round( max( sol.spos.rs( :, 1 )) - min( sol.spos.rs( :, 1 )) + sol.sz.r ) );
% sol.sample.sz.c = round( 1.00 * round( max( sol.spos.rs( :, 2 )) - min( sol.spos.rs( :, 2 )) + sol.sz.c ) );

% % round (ceiling) to nearest even
% sol.sample.sz.r = sol.sample.sz.r + mod( sol.sample.sz.r, 2 );
% sol.sample.sz.c = sol.sample.sz.c + mod( sol.sample.sz.c, 2 );
% 
% % a bit more padding in case we want to do some scan position correction
% sol.sample.sz.r = sol.sample.sz.r + 500;
% sol.sample.sz.c = sol.sample.sz.c + 500;

%}


%=================
% for zjiang202206
%=================

sol.sample.sz.r = expt.sample_sz.r;
sol.sample.sz.c = expt.sample_sz.c;

% sol.sample.sz.r = 1264;
% sol.sample.sz.c = 1264;

%=================
% for zjiang202204
%=================

% !!!!!!!! DEFINE THESE IN EXPT SETUP, ADD TO EXPT STRUCT SAMPLE PROBLEM SIZE !!!!!!!!!!!!!!

% % 180, L0339_to_L034_combined_512x512.mat
% sol.sample.sz.r = 1536;
% sol.sample.sz.c = 2688;

% % +/-5, L0342_to_L0345_combined_1024x512
% sol.sample.sz.r = 2688;
% sol.sample.sz.c = 5632;

% % +/-10, L0350_to_L0353_combined_1024x512
% sol.sample.sz.r = 2560;
% sol.sample.sz.c = 8704;

% % +/- 15, L0358_to_L0362_combined_1024x256, L0363_to_L0367_combined_1024x256
% sol.sample.sz.r = 2560;
% sol.sample.sz.c = 5760;

% % +/- 20, L0368_to_L0372_combined_1024x256, L0373_to_L0377_combined_1024x256
% sol.sample.sz.r = 2560;
% sol.sample.sz.c = 7424;

% % +/- 30, L0390_to_L0395_combined_1024x256, L0396_to_L0401_combined_1024x256
% sol.sample.sz.r = 2560;
% sol.sample.sz.c = 10240;

% % +/- 2.5, L0600_to_L0602_combined_512x512, L0603_to_L0605_combined_512x512
% % 182.5, L0606_to_L0608_combined_512x512
% % 177.5, L0609_to_L0611_combined_512x512
% sol.sample.sz.r = 1536;
% sol.sample.sz.c = 4096;

% % +/- 7.5, L0612_to_L0615_combined_1024x512, L0616_to_L0619_combined_1024x512
% % 187.5, L0620_to_L0623_combined_1024x512
% % 172.5, L0624_to_L0627_combined_1024x512
% sol.sample.sz.r = 2688;
% sol.sample.sz.c = 7040;

%==============================================================================================================================
% determine proper view region FOV ( based on probe size ) for when computing the exit wave views at a particular scan position
%==============================================================================================================================

sol.sample.vs.r = single( round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r ))));
sol.sample.vs.c = single( round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c ))));

%{

% ( r, c ) indices for each scan position:
vsr_s = transpose( round( sol.sample.vs.r - 1.0 * sol.spos.rs( :, 1 )));
vsc_s = transpose( round( sol.sample.vs.c - 1.0 * sol.spos.rs( :, 2 )));

figure; imagesc( ( vsr_s > sol.sample.sz.r ) + 2 * ( vsr_s < 1 ) )
figure; imagesc( ( vsc_s > sol.sample.sz.c ) + 2 * ( vsc_s < 1 ) )


figure; 
plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  

%}

%====================================
% book-keeping for final sample sizes
%====================================

sol.sample.sz.sz      = [ sol.sample.sz.r, sol.sample.sz.c ];
sol.sample.sz.rc      = sol.sample.sz.r * sol.sample.sz.c;
sol.sample.sz.sqrt_rc = sqrt( sol.sample.sz.rc );

%===================================
% complex valued random sample start
%===================================

sol.sample.T = rand( sol.sample.sz.r, sol.sample.sz.c ) .* exp( 0 *  1i * 2 * pi * rand( sol.sample.sz.r, sol.sample.sz.c ));
% sol.sample.T = rand( sol.sample.sz.r, sol.sample.sz.c ) + 1i * rand( sol.sample.sz.r, sol.sample.sz.c );

%================================
% Gaussian blur the sample start?
%================================

sol.sample.T = lpf_gauss( sol.sample.T, [ 0.01 * sol.sample.sz.r, 0.01 * sol.sample.sz.r ] );

%=================================
% enforce modulus and phase limits
%=================================

sol.sample.absL = 0.0 - 0.000;
sol.sample.absH = 1.0 + 0.000;
sol.sample.T   = modulus_limits_scale( sol.sample.T, [ sol.sample.absL, sol.sample.absH ] );

sol.sample.phsL = -pi * 0.01;
sol.sample.phsH = +pi * 0.01;
sol.sample.T   = phase_limits_scale( sol.sample.T, [ sol.sample.phsL, sol.sample.phsH ] );
% sol.sample.T   = phase_limits_project( sol.sample.T, [ sol.sample.phsL, sol.sample.phsH ] );

%========

% does_the_sample_array_fit_all_spos( sol, expt );

end

%====================================================================================================================================================

function does_the_sample_array_fit_all_spos( sol, expt )


    sol.spos.shifttype = 'px';
    % expt.spos.shifttype = 'subpx';


    pm = size( sol.probe.phi, 3 );

    max_P = max(max( abs(sol.probe.phi( :, :, pm ))));
    % tmp0 = single( abs(sol.probe.phi( :, :, pm )) > 0.01 * max_P );
    tmp0 = single( abs(sol.probe.phi( :, :, pm )) > 0.01 * max_P );
    % tmp9 = padarray( tmp0, double( 0.5 * [ sol.sample.sz.r - expt.sz.r, sol.sample.sz.c - expt.sz.c ] ));
    tmp9 = zeropadarray( 1 + 0 * tmp0, double( round( 0.5 * [ sol.sample.sz.r - sol.sz.r, sol.sample.sz.c - sol.sz.c ] )));

    % % just used as a visual aide to see where ~ center of probe is
    % tmp9( size( tmp9, 1 ) / 2 + 1, size( tmp9, 2 ) / 2 + 1 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 0, size( tmp9, 2 ) / 2 + 1 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 0, size( tmp9, 2 ) / 2 + 0 ) = 10;  
    % tmp9( size( tmp9, 1 ) / 2 + 1, size( tmp9, 2 ) / 2 + 0 ) = 10;  

    % tmp0 = ( abs(sol.sample.T) > 0.01 * max(abs(sol.sample.T(:))) );
    tmp55 = 0;

    for ii = 1 : sol.spos.N

    %     rj = round( -spos.rs( sio( ii ), : ));
        rj = round( +sol.spos.rs( ii, : ));

        tmp2 = circshift( round( tmp9 ), -rj );

%         [ phi, ~ ] = enforce_2DTPAsposview( 1 + 0 * sol.probe.phi( :, :, pm ), sol.sample.T, sol.sample.vs.r, sol.sample.vs.c, +rj, sol.spos.shifttype );

        tmp55 = tmp55 + tmp2 .* sol.sample.T;


        %{

        figure(666);
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1] );

    %     tmp0 = expt.meas.SI( ii ).Dneq0 .* expt.meas.SI( ii ).D + expt.meas.SI( ii ).Deq0 .* fft2( phi ) / sol.sz.sqrt_rc;
        tmp0 = fft2( phi ) / sol.sz.sqrt_rc;

    % p(1)=subplot(3,4,[1 2 5 6]);
    % p(2)=subplot(3,4,[3 4 7 8]);
    % p(3)=subplot(3,4,[9 10]);
    % p(4)=subplot(3,4,[11 12])

        subplot( 2, 3, 1 ); imagesc( log10(1 + abs( fftshift( tmp0 )))); daspect([1 1 1]); colormap jet; colorbar
        subplot( 2, 3, 4 ); imagesc( log10(1 + abs( fftshift( expt.meas.SI( ii ).D )))); daspect([1 1 1]); colormap jet; colorbar
        %subplot(122); imagesc( ( tmp1 > 0 ) .* abs( 1 + 0 * sol.sample.T ) ); daspect([1 1 1]); colormap jet
        %subplot(132); imagesc( abs(phi) > 0 ); daspect([1 1 1]); colormap jet
        subplot( 2, 3, [ 2, 5 ] ); imagescHSV( phi + 0 * tmp0 ); daspect([1 1 1]); %colormap jet
        subplot( 2, 3, [ 3, 6 ] ); imagescHSV( tmp2 .* sol.sample.T ); daspect([1 1 1]); %colormap jet
    %     tic
    %     export_fig( num2str( ii, 'test_spos_%d' ), '-r90.0' )
    %     toc
        pause( 0.01 )

        %}


    end

    figure; imagesc(abs(tmp55)>0); colormap gray; daspect([1 1 1])

    % figure; imagescHSV( ( tmp1 > 0 ) .* ( 1 + 0 * sol.sample.T )); daspect([1 1 1]); %colormap gray

    close all;


end



% 
% Z = sol.spos.rs( :, 2 ) + sol.sample.vs.c;
% sol.sample.sz.c = round( max( Z(:) ));
% 
% Z = sol.spos.rs( :, 1 ) + sol.sample.vs.r;
% sol.sample.sz.r = round( max( Z(:) ));
% 
% 
% % 
% % sol.spos.rs( :, 2 ) 

%{

close all;
figure; plot( sol.spos.rs( :, 2 ) > max( sol.sample.vs.c ) )
figure; plot( sol.spos.rs( :, 2 ) > min( sol.sample.vs.c ) )

%}

% 
% 
% 5;

