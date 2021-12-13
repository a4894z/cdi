%

%{

% codelocation = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/'; 
codelocation = '~/Documents/Science/Matlab/Code/cdi/'; 

cd( codelocation ); 
restoredefaultpath; 
addpath( genpath( pwd ));

clear; close all; testing_poisson_exwv_steplength

%}

%====================================================================================================================================================

% load /net/s8iddata/export/8-id-ECA/Analysis/atripath/sim_ptycho2DTPA_poissonsteplength.mat
load /home/ash/Documents/Science/Matlab/Code/cdi/sim_ptycho2DTPA_noisy.mat

rng( 'shuffle' )
reset( gpuDevice )

sol.spos.rs = single( sol.spos.rs );

%====================================================================================================================================================

% CHEAT CODES

a = single( 0.2 );
sol.probe.phi = ( 1 - a ) * expt.probe.phi + a * rand( [ expt.sz.sz, expt.probe.scpm.N ], 'single' );
sol.probe.scpm = expt.probe.scpm;

a = single( 0.2 );
sol.sample.T = ( 1 - a ) * expt.sample.T + a * rand( expt.sample.sz.sz, 'single' );

% sol.sample.T( abs( sol.sample.T ) > 1 ) = 1;
% sol.sample.T = lpf_gauss( sol.sample.T, [ 0.06 * expt.sample.sz.r, 0.06 * expt.sample.sz.c ] );

clear( 'a' )

%========

% sol.sample.vs.r = round( ( 0.5 * ( sol.sample.sz.r - 4 ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + 4 )));
% sol.sample.vs.c = round( ( 0.5 * ( sol.sample.sz.c - 6 ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + 6 )));

%====================================================================================================================================================

spos_N = single( round( 0.10 * sol.spos.N ));

% spos_test = sort( single( randperm( sol.spos.N, spos_N )));
spos_test = 1 : spos_N;

sample_sposview_indices = get_indices_2Dframes( sol.spos.rs( spos_test, : ), sol.sample.sz.sz, sol.sample.vs.r, sol.sample.vs.c );

ind        = uint32( sample_sposview_indices ); 
ind_offset = uint32( 0 : sol.sample.sz.rc : ( sol.sample.sz.rc * ( spos_N - 1 )));
ind_offset = ind + ind_offset;
    
%=================================================
% get sample view frames across all scan positions
%=================================================

T = reshape( sol.sample.T( ind ), [ sol.sz.sz, 1, spos_N ]);

%=========================================
% form exit wave across all scan positions
%=========================================

psi = squeeze( T .* sol.probe.phi );

%==============================================================================
% measurement constraints assuming Gaussian noise with (ignored) constant stdev
%==============================================================================

Psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;
Psi = reshape( Psi, [ sol.sz.rc, sol.probe.scpm.N, spos_N ] );

abs2_Psi = abs( Psi ) .^ 2;

% clear( 'T', 'ind', 'ind_offset', 'sample_sposview_indices' )

%====================================================================================================================================================

Nalpha = single( 40 );
alpha_test = linspace( single( 0.00 ), single( 4 ), Nalpha );

% alpha_test = single( 1 * ( 0.00 : 0.025 : 1.00 ));
% Nalpha = length( alpha_test );

%====================================================================================================================================================

use_GPU = logical( 1 );

if use_GPU == true

    expt.meas.D      = gpuArray( expt.meas.D );
    Psi              = gpuArray( Psi );
    abs2_Psi         = gpuArray( abs2_Psi );
    alpha_test       = gpuArray( alpha_test );
    expt.spos.N      = gpuArray( expt.spos.N );
    spos_N           = gpuArray( spos_N );
    expt.sz.rc       = gpuArray( expt.sz.rc );
    sol.sz.rc        = gpuArray( sol.sz.rc );
    sol.probe.scpm.N = gpuArray( sol.probe.scpm.N );
    Nalpha           = gpuArray( Nalpha );

end

%====================================================================================================================================================

I_m = reshape( expt.meas.D( :, :, spos_test ) .^ 2, [ expt.sz.rc, spos_N ] );
I_e = squeeze( sum( abs2_Psi, 2 )); % sum over scpm

xi = 1 - I_m ./ I_e;
xi( isnan( xi ) | isinf( xi ) ) = 0;

clear( 'I_e' )

%=============================================
% compute the step length by exact line search
%=============================================

tic

grad_Psi = Psi .* reshape( xi, [ sol.sz.rc, 1, spos_N ] );

abs2_Psi_alpha = abs( Psi - reshape( alpha_test, [ 1, 1, 1, length( alpha_test ) ] ) .* grad_Psi ) .^ 2;

A_i_q_alpha_p = abs2_Psi_alpha + sum( abs2_Psi, 2 ) - abs2_Psi;

poisson_cost = squeeze( sum( A_i_q_alpha_p - reshape( I_m, [ expt.sz.rc, 1, spos_N, 1 ] ) .* log( A_i_q_alpha_p ), 1 ));

toc

%========

poisson_cost = permute( poisson_cost, [ 2, 3, 1 ] );
poisson_cost = poisson_cost - min( poisson_cost, [], 2 );

kk = 1;

for pp = 1 : sol.probe.scpm.N

    figure( 666 );
    subplot( 1, 3, kk )
%     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
    imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + abs( poisson_cost( :, :, pp ) )))
    title( num2str( [ pp, sol.probe.scpm.occ( pp ) ], 'pp = %d, occupancy = %0.6f' ))
    grid on
    colormap turbo
    xlabel('\alpha')
    ylabel('spos index')
    
    kk = kk + 1;

end

%========

clear( 'pp', 'kk', 'A_i_q_alpha_p', 'abs2_Psi_alpha' )
 
%================================================
% compute the step length by 1st order optimality
%================================================

xi_abs_Psi2 = reshape( xi, [ sol.sz.rc, 1, spos_N ] ) .* abs2_Psi;

tic

xi_alpha_minus_one = xi .* reshape( alpha_test, [ 1, 1, Nalpha ] ) - 1;
lhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* reshape( xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] ), 1 ));

numer = reshape( I_m .* xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] );
denom = abs2_Psi .* reshape( xi_alpha_minus_one .^ 2, [ sol.sz.rc, 1, spos_N, Nalpha ] ) + sum( abs2_Psi, 2 ) - abs2_Psi ;
rhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));

f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );

toc

%========

f_eq_0 = permute( f_eq_0, [ 2, 1, 3 ] );
f_eq_0 = f_eq_0 - min( f_eq_0, [], 3 );

kk = 1;

for pp = 1 : sol.probe.scpm.N
    
    figure( 667 );
    subplot( 1, 3, kk )
%     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
    imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0( :, pp, : )) ))
%     imagesc( log10( 1 + abs( squeeze( f_eq_0( :, pp, : )) )))
    title( num2str( [ pp, sol.probe.scpm.occ( pp ) ], 'pp = %d, occupancy = %0.6f' ))
    grid on
    colormap turbo
    xlabel('\alpha')
    ylabel('spos index')
    
    kk = kk + 1;
    
end

% [ ~, II ] = min( f_eq_0, [], 3 );
% 
% alpha_scpm_vs_spos = alpha_test( II );

%============================================
% actually perform gradient descent iteration
%============================================

[ ~, II ] = min( f_eq_0, [], 3 );
alpha_opt = transpose( alpha_test( II ));
alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ] );

Psi_opt   = reshape( Psi - alpha_opt .* grad_Psi, [ sol.sz.sz, sol.probe.scpm.N, spos_N ] );

psi_opt   = fftshift( fftshift( ifft( ifft( Psi_opt, [], 1 ), [], 2 ), 1 ), 2 ) * sol.sz.sqrt_rc;

% figure; imagesc( abs( psi(:,:,3,44) ), [0, 100 ]); 
% figure; imagesc( abs( psi_opt(:,:,3,44) ), [0, 100 ])

% %=================================================================
% % sample/probe update direction scaling for approximation in Eq 21
% %=================================================================
% 
% chi = psi_opt - psi;
% 
% Delta_P = chi .* conj( T );
% Delta_O = chi .* conj( sol.probe.phi );
% 
% % t1 = T .* Delta_P;
% % t2 = sol.probe.phi .* Delta_O;
% % t3 = Delta_P .* Delta_O;
% %
% % figure; imagesc( abs( chi( :, :, 3, 44 ) )); colorbar; title('chi');
% % figure; imagesc( abs( psi( :, :, 3, 44 ) )); colorbar; title('psi');
% % figure; imagesc( abs( t1( :, :, 3, 44 ) ));  colorbar; title('T .* Delta\_P')
% % figure; imagesc( abs( t2( :, :, 3, 44 ) ));  colorbar; title('sol.probe.phi .* Delta\_O')
% % figure; imagesc( abs( t3( :, :, 3, 44 ) ));  colorbar; title('Delta\_P .* Delta\_O')
% 
% figure; imagesc( abs(Delta_P(:,:,3,ss)) );    daspect([1 1 1]); colormap turbo; colorbar; title('Delta\_P')
% figure; imagesc( abs(Delta_O(:,:,3,ss)) );    daspect([1 1 1]); colormap turbo; colorbar; title('Delta\_O')
% figure; imagesc( abs(chi(:,:,3,ss)) );        daspect([1 1 1]); colormap turbo; colorbar; title('chi')
% figure; imagesc( abs(sol.probe.phi(:,:,3)) ); daspect([1 1 1]); colormap turbo; colorbar; title('P')
% figure; imagesc( abs(T(:,:, 1, ss )) );       daspect([1 1 1]); colormap turbo; colorbar; title('O')
% 
% %========
%         
% pp = 1;
% ss = 23;
% 
% t1 = chi( :, :, pp, ss );
% t2 = sol.probe.phi( :, :, pp ) .* Delta_O( :, :, pp, ss );
% t3 = T( :, :, :, ss )          .* Delta_P( :, :, pp, ss );
% t4 = Delta_P( :, :, pp, ss )   .* Delta_O( :, :, pp, ss );
% 
% alpha_P_test = linspace( -5.0, 5.0, 100 );
% alpha_O_test = linspace( -0.05, 0.05, 100 );     % for pp = 1
% % alpha_O_test = linspace( -0.05, 0.05, 100 );       % for pp = 2
% % alpha_O_test = linspace( -0.0005, 0.0005, 100 );   % for pp = 3
% 
% L_alphaP_alphaO        = zeros( length( alpha_P_test ), length( alpha_O_test ), 'single' );
% L_alphaP_alphaO_approx = zeros( length( alpha_P_test ), length( alpha_O_test ), 'single' );
% % cost_dropped_term      = zeros( length( alpha_P_test ), length( alpha_O_test ), 'single' );
% 
% for PP = 1 : length( alpha_P_test )
%     
%     for OO = 1 : length( alpha_O_test )
%         
%         tmp0 = t1 - alpha_O_test( OO ) * t2 ...
%                   - alpha_P_test( PP ) * t3;
%         
%         tmp1 = -alpha_O_test( OO ) * alpha_P_test( PP ) * t4;
%         
%         tmp2 = tmp0 + tmp1;
%         
%         L_alphaP_alphaO_approx( PP, OO ) = sum( abs( tmp0( : )) .^ 2 );
% %         cost_dropped_term( PP, OO )      = sum( abs( tmp1( : )) .^ 2 );
%         L_alphaP_alphaO( PP, OO )        = sum( abs( tmp2( : )) .^ 2 );
% 
%         
%     end
%         
% end
% 
% figure; 
% imagesc( alpha_O_test, alpha_P_test, log10( L_alphaP_alphaO ))
% set( gca, 'ydir', 'normal') 
% grid on
% colormap turbo
% xlabel('alpha\_O')
% ylabel('alpha\_P')
% title('Least Squares Cost Function')
% 
% figure; 
% imagesc( alpha_O_test, alpha_P_test, log10( L_alphaP_alphaO_approx ))
% set( gca, 'ydir', 'normal') 
% grid on
% colormap turbo
% xlabel('alpha\_O')
% ylabel('alpha\_P')
% title('Approximate Least Squares Cost Function')
% 
% figure; 
% imagesc( alpha_O_test, alpha_P_test, abs( log10( L_alphaP_alphaO_approx ) - log10( L_alphaP_alphaO )) )
% set( gca, 'ydir', 'normal') 
% grid on
% colormap turbo
% xlabel('alpha\_O')
% ylabel('alpha\_P')
% title('Difference')
% 
% return


%===========================================================================
% using results from ELS, determine if nonlinear eqn is positive or negative
%                       if +, then decrease step size
%                       if -, then increase step size
%===========================================================================

% give updated exitwave a small perturbation to kick it out of current local solution
psi = psi_opt + 0.1 * rand( size( psi_opt )); 

%========

Psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;
Psi = reshape( Psi, [ sol.sz.rc, sol.probe.scpm.N, spos_N ] );

abs2_Psi = abs( Psi ) .^ 2;

I_e = squeeze( sum( abs2_Psi, 2 )); % sum over scpm

xi = 1 - I_m ./ I_e;
xi( isnan( xi ) | isinf( xi ) ) = 0;

clear( 'I_e' )

%========

###
xi_abs_Psi2 = reshape( xi, [ sol.sz.rc, 1, spos_N ] ) .* abs2_Psi;

xi_alpha_minus_one = xi .* reshape( alpha_test, [ 1, 1, Nalpha ] ) - 1;

lhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* reshape( xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] ), 1 ));

numer              = reshape( I_m .* xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] );
denom              = abs2_Psi .* reshape( xi_alpha_minus_one .^ 2, [ sol.sz.rc, 1, spos_N, Nalpha ] ) + sum( abs2_Psi, 2 ) - abs2_Psi ;
rhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));

f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );











return


%======================================================
% increasing step until "zero" is found, using for loop
%======================================================

xi  = reshape( xi,  [ sol.sz.rc, 1, spos_N ] );
I_m = reshape( I_m, [ sol.sz.rc, 1, spos_N ] );

xi_abs_Psi2 = xi .* abs2_Psi;

sum_pprime_neq_p = ( sum( abs2_Psi, 2 ) - abs2_Psi );


% alpha_test = linspace( 0.0, 2, 10 );
alpha_test_2 = alpha_test;
spos_test_2  = spos_test;

alpha_found = zeros( spos_N, 1, 'single' );

%========

% test a step length to see if the step length equation is positive

pp = 3;

tic

% for aa = 1 : 5
for aa = 1 : length( alpha_test_2 )
    
    %========
    
    xi_alpha_minus_one = xi .* alpha_test_2( aa ) - 1;

    lhs_steplength_eqn = squeeze( sum( xi_abs_Psi2( :, pp, : ) .* xi_alpha_minus_one, 1 ));

    numer = xi_alpha_minus_one .* I_m;
    denom = abs2_Psi( :, pp, : ) .* ( xi_alpha_minus_one .^ 2 ) + sum_pprime_neq_p( :, pp, : );
    rhs_steplength_eqn = squeeze( sum( xi_abs_Psi2( :, pp, : ) .* ( numer ./ denom ), 1 ));
    
    %========
    
    goofy_nonlinear_eqn_equal_0 = lhs_steplength_eqn - rhs_steplength_eqn;

    found_a_step         = ( goofy_nonlinear_eqn_equal_0 > 0 );
    step_found_spos_indx = spos_test_2( found_a_step );

    alpha_found( step_found_spos_indx ) =  alpha_test_2( aa );
        
%     if isempty( step_found_spos_indx ), continue; end

    %========

    % remove scan positions from the process since we've found a step length for them
    
    spos_test_2( found_a_step ) = [];
    if isempty( spos_test_2 ), break; end
    
    xi(               :, :, found_a_step ) = [];
    I_m(              :, :, found_a_step ) = [];
    xi_abs_Psi2(      :, :, found_a_step ) = [];
    sum_pprime_neq_p( :, :, found_a_step ) = [];
    abs2_Psi(         :, :, found_a_step ) = [];

end

toc

tmp0 = squeeze( f_eq_0( :, pp, : ));
[ ~, II ] = min( tmp0, [], 2 );

max_alpha = max( alpha_test );
max_alpha = 1.2 * gather( max_alpha );

figure; 
hold on
plot( alpha_test( II ), 'linewidth', 3, 'color', [ 0.0, 0.0, 0.0 ] )
plot( alpha_found, '-.', 'linewidth', 2, 'color', [ 0.0, 0.7, 0.0 ] )
hold off
ylim( [ 0.0, max_alpha ] )
grid on
xlabel('scan position index')
ylabel('step length for exitwave')
title( num2str( [ pp, sol.probe.scpm.occ( pp ) ], 'pp = %d, occupancy = %0.6f' ))

return

%====================================================================================================================================================

