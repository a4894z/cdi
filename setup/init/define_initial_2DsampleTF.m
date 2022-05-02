function [ sol, expt ] = define_initial_2DsampleTF( sol, expt )

%=====================================
% load previous sample here if desired
%=====================================

% [ sol, expt ] = load_previous_2DsampleTF( sol, expt );
% return

%==========================
% array size for the sample
%==========================

sol.sample.sz.r = round( 1.00 * round( max( sol.spos.rs( :, 1 )) - min( sol.spos.rs( :, 1 )) + sol.sz.r ) );
sol.sample.sz.c = round( 1.00 * round( max( sol.spos.rs( :, 2 )) - min( sol.spos.rs( :, 2 )) + sol.sz.c ) );

%=======================================================================
% a bit more padding in case we want to do some scan position correction
%=======================================================================

sol.sample.sz.r = sol.sample.sz.r + 20;
sol.sample.sz.c = sol.sample.sz.c + 40;

%================================
% round (ceiling) to nearest even
%================================

sol.sample.sz.r = sol.sample.sz.r + mod( sol.sample.sz.r, 2 );
sol.sample.sz.c = sol.sample.sz.c + mod( sol.sample.sz.c, 2 );

%==============================================================================================================================
% determine proper view region FOV ( based on probe size ) for when computing the exit wave views at a particular scan position
%==============================================================================================================================

sol.sample.vs.r = single( round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r ))));
sol.sample.vs.c = single( round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c ))));

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

%========

% does_the_sample_array_fit_all_spos( sol, expt );

%================================
% Gaussian blur the sample start?
%================================

sol.sample.T = lpf_gauss( sol.sample.T, [ 0.005 * sol.sample.sz.r, 0.005 * sol.sample.sz.r ] );

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

