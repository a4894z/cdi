%

%{

codelocation = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/'; 
% codelocation = '~/Documents/Science/Matlab/Code/cdi/'; 

cd( codelocation ); 
restoredefaultpath; 
addpath( genpath( pwd ));

clear; close all; testing_poisson_exwv_steplength



%}



%====================================================================================================================================================

load /net/s8iddata/export/8-id-ECA/Analysis/atripath/sim_ptycho2DTPA_poissonsteplength.mat
% load /home/ash/Documents/Science/Matlab/Code/cdi/sim_ptycho2DTPA.mat

rng( 'shuffle' )
reset( gpuDevice )

sol.spos.rs = single( sol.spos.rs );

%====================================================================================================================================================

% CHEAT CODES

a = single( 0.2 );
sol.probe.phi = ( 1 - a ) * expt.probe.phi + a * rand( [ expt.sz.sz, expt.probe.scpm.N ], 'single' );
sol.probe.scpm = expt.probe.scpm;

a = single( 0.4 );
sol.sample.T = ( 1 - a ) * expt.sample.T + a * rand( expt.sample.sz.sz, 'single' );

% sol.sample.T( abs( sol.sample.T ) > 1 ) = 1;
% sol.sample.T = lpf_gauss( sol.sample.T, [ 0.06 * expt.sample.sz.r, 0.06 * expt.sample.sz.c ] );

clear( 'a' )

%========

% sol.sample.vs.r = round( ( 0.5 * ( sol.sample.sz.r - 4 ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + 4 )));
% sol.sample.vs.c = round( ( 0.5 * ( sol.sample.sz.c - 6 ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + 6 )));

%====================================================================================================================================================

spos_N = single( round( 0.25 * sol.spos.N ));

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

psi = T .* sol.probe.phi;

psi = squeeze( psi );

%==============================================================================
% measurement constraints assuming Gaussian noise with (ignored) constant stdev
%==============================================================================

Psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sol.sz.sqrt_rc;
Psi = reshape( Psi, [ sol.sz.rc, sol.probe.scpm.N, spos_N ] );

abs2_Psi = abs( Psi ) .^ 2;

clear( 'psi', 'T', 'ind', 'ind_offset', 'sample_sposview_indices' )

%====================================================================================================================================================

Nalpha = single( 50 );

alpha_test = linspace( single( 0.00 ), single( 3 ), Nalpha );

%====================================================================================================================================================

use_GPU = logical( 1 );

if use_GPU == true

    expt.meas.D      = gpuArray( expt.meas.D );
    Psi              = gpuArray( Psi );
    alpha_test       = gpuArray( alpha_test );
    expt.spos.N      = gpuArray( expt.spos.N );
    spos_N           = gpuArray( spos_N );
    expt.sz.rc       = gpuArray( expt.sz.rc );
    sol.sz.rc        = gpuArray( sol.sz.rc );
    sol.probe.scpm.N = gpuArray( sol.probe.scpm.N );
    Nalpha           = gpuArray( Nalpha );

end

%====================================================================================================================================================

I_m      = reshape( expt.meas.D( :, :, spos_test ) .^ 2, [ expt.sz.rc, spos_N ] );

abs_Psi2 = abs( Psi ) .^ 2;
I_e      = squeeze( sum( abs_Psi2, 2 )); % sum over scpm

xi = 1 - I_m ./ I_e;
xi( isnan( xi ) | isinf( xi ) ) = 0;

% clear( 'I_e' )

%=============================================
% compute the step length by exact line search
%=============================================

tic

grad_Psi = Psi .* reshape( xi, [ sol.sz.rc, 1, spos_N ] );

% abs2_Psi_alpha = gpuArray.zeros( [ sol.sz.rc, sol.probe.scpm.N, spos_N, Nalpha ], 'single' );
% 
% for aa =  1 : Nalpha
%     
%     abs2_Psi_alpha( :, :, :, aa ) = abs( Psi - alpha_test( aa ) .* grad_Psi ) .^ 2;
%     % I_e_alpha = abs( Psi - alpha_test .* grad_Psi ) .^ 2;
% 
% end
% 
% clear( 'aa' )

abs2_Psi_alpha = abs( Psi - reshape( alpha_test, [ 1, 1, 1, length( alpha_test ) ] ) .* grad_Psi ) .^ 2;

%========

% tic
% 
% poisson_cost = gpuArray.zeros( spos_N, Nalpha, sol.probe.scpm.N, 'single' );
% 
% kk = 1;
% 
% for pp = 1 : sol.probe.scpm.N
% 
%     A_i_q_alpha_p = squeeze( abs2_Psi_alpha( :, pp, :, : ) + sum( abs2_Psi, 2 ) - abs2_Psi( :, pp, : ) );
% 
%     poisson_cost( :, :, pp ) = squeeze( sum( A_i_q_alpha_p - I_m .* log( A_i_q_alpha_p ), 1 ));     % sum over q
% 
% %     poisson_cost( :, :, pp ) = poisson_cost( :, :, pp ) - min( poisson_cost( :, :, pp ), [], 2 );
% % 
% %     figure( 666 );
% %     subplot( 1, 3, kk )
% % %     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
% %     imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + abs( poisson_cost( :, :, pp ) )))
% %     title( num2str( pp, 'pp = %d' ))
% %     
% %     kk = kk + 1;
% 
% end
% 
% toc


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

% clear( 'pp', 'kk', 'A_i_q_alpha_p', 'grad_Psi', 'abs2_Psi_alpha' )
 
%================================================
% compute the step length by 1st order optimality
%================================================

xi_abs_Psi2 = reshape( xi, [ sol.sz.rc, 1, spos_N ] ) .* abs_Psi2;

% tic
% 
% lhs_steplength_eqn_B = gpuArray.zeros( spos_N, sol.probe.scpm.N, Nalpha, 'single' );
% rhs_steplength_eqn_B = gpuArray.zeros( spos_N, sol.probe.scpm.N, Nalpha, 'single' );
% 
% % xi          = reshape( xi, [ sol.sz.rc, 1, spos_N ] );
% % xi_abs_Psi2 = xi .* abs_Psi2;
% 
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     for aa =  1 : Nalpha
%         
%         tmp0 = ( alpha_test( aa ) * xi - 1 );
%         tmp2 = squeeze( xi_abs_Psi2( :, pp, : ));
%         
%         lhs_steplength_eqn_B( :, pp, aa ) = sum( tmp2 .* tmp0, 1 );
%         
%         tmp1 = squeeze( abs_Psi2( :, pp, : ));
%         
%         numer = I_m .* tmp0;
%         denom = tmp1 .* ( tmp0 .^ 2 ) + squeeze( sum( abs2_Psi, 2 ) - abs2_Psi( :, pp, : ));
%         
%         rhs_steplength_eqn_B( :, pp, aa ) = sum( tmp2 .* ( numer ./ denom ), 1 );
%     
%     end
% 
% end
% 
% toc

tic

xi_alpha_minus_one = xi .* reshape( alpha_test, [ 1, 1, Nalpha ] ) - 1;
lhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* reshape( xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] ), 1 ));

numer = reshape( I_m .* xi_alpha_minus_one, [ sol.sz.rc, 1, spos_N, Nalpha ] );
denom = abs_Psi2 .* reshape( xi_alpha_minus_one .^ 2, [ sol.sz.rc, 1, spos_N, Nalpha ] ) + sum( abs2_Psi, 2 ) - abs2_Psi ;
rhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));

% norm( lhs_steplength_eqn_B(:) - lhs_steplength_eqn(:) )
% norm( rhs_steplength_eqn_B(:) - rhs_steplength_eqn(:) )

f_eq_0   = abs( lhs_steplength_eqn - rhs_steplength_eqn   );
% f_eq_0_B = abs( lhs_steplength_eqn_B - rhs_steplength_eqn_B );

toc

f_eq_0 = permute( f_eq_0, [ 2, 1, 3 ] );
f_eq_0 = f_eq_0 - min( f_eq_0, [], 3 );
% f_eq_0_B = f_eq_0_B - min( f_eq_0_B, [], 3 );

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

% kk = 1;
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     figure( 668 );
%     subplot( 1, 3, kk )
%     imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0_B( :, pp, : )) ))
%     title( num2str( pp, 'pp = %d' ))
%     grid on
%     colormap turbo
%     
%     kk = kk + 1;
%     
% end


% kk = 1;
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     figure( 669 );
%     subplot( 1, 3, kk )
% %     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
%     imagesc( alpha_test, 1 : length( spos_test ), squeeze( f_eq_0_B( :, pp, : )) - squeeze( f_eq_0( :, pp, : )) )
% %     imagesc( log10( 1 + abs( squeeze( f_eq_0( :, pp, : )) )))
%     title( num2str( pp, 'pp = %d' ))
%     grid on
%     colormap turbo
%     
%     kk = kk + 1;
%     
% end

%============================================
% actually perform gradient descent iteration
%============================================

% Psi - reshape( alpha_test, [ 1, 1, 1, length( alpha_test ) ] ) .* grad_Psi


%====================================================================================================================================================













% 
% % tic
% 
% I_m      = reshape( expt.meas.D( :, :, spos_test ) .^ 2, [ expt.sz.rc, spos_N ] );
% 
% abs_Psi2 = abs( Psi ) .^ 2;
% I_e      = squeeze( sum( abs_Psi2, 2 )); % sum over scpm
% 
% xi = 1 - I_m ./ I_e;
% xi( isnan( xi ) | isinf( xi ) ) = 0;
% 
% %=============================================
% % compute the step length by exact line search
% %=============================================
% 
% tic
% 
% grad_Psi = Psi .* reshape( xi, [ sol.sz.rc, 1, spos_N ] );
% 
% abs2_Psi_alpha = gpuArray.zeros( [ sol.sz.rc, sol.probe.scpm.N, spos_N, Nalpha ], 'single' );
% 
% for aa =  1 : Nalpha
%     
%     abs2_Psi_alpha( :, :, :, aa ) = abs( Psi - alpha_test( aa ) .* grad_Psi ) .^ 2;
%     % I_e_alpha = abs( Psi - alpha_test .* grad_Psi ) .^ 2;
% 
% end
% 
% clear( 'aa' )
% 
% %========
% 
% poisson_cost = gpuArray.zeros( spos_N, Nalpha, sol.probe.scpm.N, 'single' );
% 
% kk = 1;
% 
% for pp = 1 : sol.probe.scpm.N
% 
%     A_i_q_alpha_p = squeeze( abs2_Psi_alpha( :, pp, :, : ) + sum( abs2_Psi, 2 ) - abs2_Psi( :, pp, : ) );
% 
%     poisson_cost( :, :, pp ) = squeeze( sum( A_i_q_alpha_p - I_m .* log( A_i_q_alpha_p ), 1 ));     % sum over q
% 
% %     poisson_cost( :, :, pp ) = poisson_cost( :, :, pp ) - min( poisson_cost( :, :, pp ), [], 2 );
% % 
% %     figure( 666 );
% %     subplot( 1, 3, kk )
% % %     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
% %     imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + abs( poisson_cost( :, :, pp ) )))
% %     title( num2str( pp, 'pp = %d' ))
% %     
% %     kk = kk + 1;
% 
% end
% 
% toc
% 
% %========
% 
% poisson_cost = poisson_cost - min( poisson_cost, [], 2 );
% 
% for pp = 1 : sol.probe.scpm.N
% 
%     figure( 666 );
%     subplot( 1, 3, kk )
% %     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
%     imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + abs( poisson_cost( :, :, pp ) )))
%     title( num2str( pp, 'pp = %d' ))
%     grid on
%     colormap turbo
%     
%     kk = kk + 1;
% 
% end
% 
% %========
% 
% clear( 'pp', 'kk', 'A_i_q_alpha_p', 'grad_Psi', 'abs2_Psi_alpha' )
%  
% %================================================
% % compute the step length by 1st order optimality
% %================================================
% 
% tic
% 
% lhs_steplength_eqn = gpuArray.zeros( spos_N, sol.probe.scpm.N, Nalpha, 'single' );
% rhs_steplength_eqn = gpuArray.zeros( spos_N, sol.probe.scpm.N, Nalpha, 'single' );
% 
% % xi          = reshape( xi, [ sol.sz.rc, 1, spos_N ] );
% % xi_abs_Psi2 = xi .* abs_Psi2;
% 
% xi_abs_Psi2 = reshape( xi, [ sol.sz.rc, 1, spos_N ] ) .* abs_Psi2;
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     for aa =  1 : Nalpha
%         
%         tmp0 = ( alpha_test( aa ) * xi - 1 );
%         tmp2 = squeeze( xi_abs_Psi2( :, pp, : ));
%         
%         lhs_steplength_eqn( :, pp, aa ) = sum( tmp2 .* tmp0, 1 );
%         
%         tmp1 = squeeze( abs_Psi2( :, pp, : ));
%         
%         numer = I_m .* tmp0;
%         denom = tmp1 .* ( tmp0 .^ 2 ) + squeeze( sum( abs2_Psi, 2 ) - abs2_Psi( :, pp, : ));
%         
%         rhs_steplength_eqn( :, pp, aa ) = sum( tmp2 .* ( numer ./ denom ), 1 );
%     
%     end
%     
% end
% 
% f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );
% 
% toc
% 
% f_eq_0 = f_eq_0 - min( f_eq_0, [], 3 );
% 
% 
% kk = 1;
% 
% for pp = 1 : sol.probe.scpm.N
%     
%     figure( 667 );
%     subplot( 1, 3, kk )
% %     imagesc( alpha_test, spos_test, log10( 1 + abs( poisson_cost )))
%     imagesc( alpha_test, 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0( :, pp, : )) ))
% %     imagesc( log10( 1 + abs( squeeze( f_eq_0( :, pp, : )) )))
%     title( num2str( pp, 'pp = %d' ))
%     grid on
%     colormap turbo
%     
%     kk = kk + 1;
%     
% end

%====================================================================================================================================================



