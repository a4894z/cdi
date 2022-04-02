function [ psi_opt, alpha_opt ] = exitwave_update_2DTPA_poisson_v2( phi,       ...
                                                                 T0,         ...
                                                                 ind,        ...
                                                                 sz,         ...
                                                                 Nspos,      ...
                                                                 Nscpm,      ...
                                                                 alpha_prev, ...
                                                                 rc,         ... 
                                                                 sqrt_rc,    ...
                                                                 I_m,        ...
                                                                 meas_Deq0,  ...
                                                                 measLPF )

%========================================
% get sample frames at each scan position
%========================================

T = reshape( T0( ind ), [ sz, 1, Nspos ]);

%===========================
% form exitwaves using 2DTPA
%===========================

psi = T .* phi;

%=====================================
% propagate exitwave modes to detector
%=====================================

Psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / sqrt_rc;
Psi = reshape( Psi, [ rc, Nscpm, Nspos ] );
Psi = permute( Psi, [ 1, 3, 2 ] );

%==================================================================
% intensity for current sample/probe solutions at measurement plane
%==================================================================

abs2_Psi = abs( Psi ) .^ 2;

I_e = squeeze( sum( abs2_Psi, 3 ));        % sum over scpm

%==================================================
% intensity for measurements over current minibatch
%==================================================

xi = 1 - I_m ./ I_e;

% xi( isnan( xi ) | isinf( xi ) ) = 0;

chknan = isnan( xi );
chkinf = isinf( xi );

if any( chknan( : ) ), xi( chknan ) = 0; end
if any( chkinf( : ) ), xi( chkinf ) = 0; end






% xi( abs( xi ) > 2 ) = 0;

%{

close all

tmp2 = reshape( I_e, [ sz, Nspos ] );  
tmp3 = reshape( I_m, [ sz, Nspos ] );  
tmp5 = reshape( xi, [ sz, Nspos ] );  

figure; imagesc( log10( 1 + abs( fftshift( tmp2(:,:,44) )))); title( 'I\_e' )
figure; imagesc( log10( 1 + abs( fftshift( tmp3(:,:,44) )))); title( 'I\_m' )
figure; imagesc( log10( 1 + abs( fftshift( tmp5(:,:,44) )))); title( 'xi' )

%}

%================================================
% compute the step length by 1st order optimality
%================================================

% xi_abs_Psi2        = xi .* abs2_Psi;
% xi_alpha_minus_one = xi .* alpha_prev - 1;
% 
% %========
% 
% lhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* xi_alpha_minus_one, 1 ));    % sum over image pixels
% 
% %========
% 
% numer = I_m .* xi_alpha_minus_one;
% denom = abs2_Psi .* xi_alpha_minus_one .^ 2 + I_e - abs2_Psi;
% 
% rhs_steplength_eqn = squeeze( sum( xi_abs_Psi2 .* ( numer ./ denom ), 1 ));    % sum over image pixels
% 
% %========
% 
% f_eq_0 = abs( lhs_steplength_eqn - rhs_steplength_eqn );
% 
% 
% 
% 
% 
% 
% 
% 
% %{
% 
% kk = 1;
% 
% for pp = 1 : Nscpm
% 
%     figure( 667 );
%     subplot( 1, double(gather( Nscpm )), kk )
%     
%     imagesc( 1 : 3, 1 : Nspos, squeeze( log10( 1 + f_eq_0( :, pp, : ) )))
% %     imagesc( squeeze( I_e( 1, 1, pp, : )), 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0n( :, pp, : )) ))
% %     imagesc( I_e( :, pp ), 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0n( :, pp, : )) ))
%     
% 
% 
% %     
% %     hold on
% %     plot( alpha_opt( 1, :, pp ), 1 : Nspos, 'r.', 'MarkerSize', 5 );
% % %     plot( alpha_opt, 1 : size( f_eq_0, 1 ), 'r.', 'MarkerSize', 5 );
% 
%     hold off
%     
%     grid on
%     colormap bone
%     xlabel('Test \alpha Index')
%     ylabel('spos index')
%     
%     kk = kk + 1;
%     
%     drawnow
% 
% end
% 
% %}
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% [ ~, II ] = min( f_eq_0, [], 3 );
% % [ ~, II ] = min( f_eq_0, [], 2 );
% 
% alpha_opt = gpuArray.zeros( Nspos, Nscpm, 'single' );
% 
% 
% 
% 
% 
% for ss = 1 : Nspos
%     
%     for pp = 1 : Nscpm
% 
%         alpha_opt( ss, pp ) = alpha_prev( 1, ss, pp, II( ss, pp ) );
% 
%     end
% 
% end
% 
% 
% 
% 
% 
% % for pp = 1 : Nscpm
% % 
% %     alpha_opt( :, pp ) = alpha_prev( 1, :, pp, II( :, pp ) );
% % 
% % end
% 
% 
% 
% 
% alpha_opt = reshape( alpha_opt, [ 1, size( alpha_opt ) ] );






alpha_opt = alpha_prev;





%{






kk = 1;

for pp = 1 : Nscpm

    figure( 667 );
    subplot( 1, double(gather( Nscpm )), kk )
    
    imagesc( 1 : 3, 1 : Nspos, squeeze( log10( 1 + f_eq_0( :, pp, : ) )))
%     imagesc( squeeze( I_e( 1, 1, pp, : )), 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0n( :, pp, : )) ))
%     imagesc( I_e( :, pp ), 1 : length( spos_test ), log10( 1 + squeeze( f_eq_0n( :, pp, : )) ))
    

    tmp0 = 1 : 3;
    
    hold on
    plot( II( :, pp ), 1 : Nspos, 'r.', 'MarkerSize', 5 );
%     plot( alpha_opt, 1 : size( f_eq_0, 1 ), 'r.', 'MarkerSize', 5 );

    hold off
    
    grid on
    colormap bone
    xlabel('\alpha Index')
    ylabel('spos index')
    
    kk = kk + 1;
    
    drawnow

end

%}







%====================================================================================================================================================

Psi_opt = reshape( Psi .* ( 1 - alpha_opt .* xi ), [ sz, Nspos, Nscpm ] );



% grad_Psi = Psi .* xi;
% Psi_opt   = reshape( Psi - alpha_opt .* grad_Psi, [ sz, Nscpm, Nspos ] );
% Psi_opt   = reshape( Psi - 0.25 .* grad_Psi, [ sz, Nscpm, Nspos ] );

psi_opt   = fftshift( fftshift( ifft( ifft( Psi_opt, [], 1 ), [], 2 ), 1 ), 2 ) * sqrt_rc;

psi_opt = permute( psi_opt, [ 1, 2, 4, 3 ]);

5;

%{

close all

figure; imagesc( abs(psi_opt(:,:,1,44)  ) ) 
figure; imagesc( abs(psi_opt(:,:,2,44)  ) ) 
figure; imagesc( abs(psi_opt(:,:,3,44)  ) ) 

figure; imagesc( abs(psi(:,:,1,44)  ) ) 
figure; imagesc( abs(psi(:,:,2,44)  ) ) 
figure; imagesc( abs(psi(:,:,3,44)  ) ) 







tmp0 = reshape( Psi, [ sz, Nscpm, Nspos ] );  
tmp1 = reshape( grad_Psi, [ sz, Nscpm, Nspos ] );
tmp2 = reshape( I_e, [ sz, Nspos ] );  
tmp3 = reshape( I_m, [ sz, Nspos ] );  
tmp4 = reshape( meas_Deq0, [ sz, Nspos ] );  

close all

figure; imagesc( log10( 1 + abs( fftshift( tmp0(:,:,44) )))); title( 'Psi' )
figure; imagesc( log10( 1 + abs( fftshift( tmp1(:,:,44) )))); title( 'grad\_Psi' )
figure; imagesc( log10( 1 + abs( fftshift( tmp2(:,:,44) )))); title( 'I\_e' )
figure; imagesc( log10( 1 + abs( fftshift( tmp3(:,:,44) )))); title( 'I\_m' )
figure; imagesc( log10( 1 + abs( fftshift( tmp4(:,:,44) )))); title( 'meas\_Deq0' )

figure; imagesc( log10( 1 + abs( fftshift( not(tmp4(:,:,44)) .* tmp1(:,:,44) )))); 


tmp5 = reshape( xi, [ sz, Nspos ] );  
figure; imagesc( fftshift(   tmp5(:,:,44) )); title( 'xi' )




tmp6 = reshape( 1 - alpha_opt .* xi, [ sz, Nscpm, Nspos ] );
figure; imagesc( log10( 1 + abs( fftshift(   tmp6(:,:,3, 44) )))); title( 'exitwave scaling update' )






%}

