function V = enforce_2DTPAmeas( phi, meas, meas_Deq0, measLPF, sol )

% sol.algos.measupdate = { 'std', 'gauss', 'poiss', 'ring' };



% sol.scpm.N             % size( phi, 3 )
% sol.sz.r               % size( phi, 1 )   
% sol.sz.sqrt_rc         % sqrt_rc = sqrt( size( phi, 1 ) * size( phi, 2 ))
% meas.SI
% meas.SIeq0

%==================================================================================================

% the measurement projection operator is
% 
% r is the spatial position coordinate
% q is the spatial frequency coordinate
% k is index for probe modes
% l is index for sample modes
% s is scan position index
% 
% psi_j( r, k, l ) = ifft2[ sqrt(I_j( q ))  * phi_tilde_j( r, k, l ) / sqrt( sum_{k, l}[ | phi_tilde_j( r, k, l ) |^2 ] )  ]

%==================================================================================================

% if isempty( missing_data )
%     if ~exist( missing_data, 'in' )
%         missing_data.in = ( measurement == 0 ); 
%         missing_data.out = not( missing_data.in ); 
%         warning('Assigning missing data region in measurement projection operator on the fly, very inefficient');
%     end
% end

% try
%    missing_data;
% catch ME
%    missing_data.in = ( measurement == 0 ); 
%    missing_data.out = not( missing_data.in ); 
%    warning('Assigning missing data region in measurement projection operator on the fly, very inefficient');
% end

%==================================================================================================
% calculate the current iterate for exit wave field intensity 
% at the measurement plane z3 over all probe modes:
%==================================================================================================

% form term sum_p | phi_tilde_{ s, p } |^2
% TODO: add option for fresnel propagator here

% % init storage
% V = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );
% tmp0 = zeros( sol.sz.r, sol.sz.c, 'single' );
% 
% for pp = 1 : sol.probe.scpm.N             % size( phi, 3 )
% 
%     V( :, :, pp ) = fft2( fftshift( phi( :, :, pp ))) / sol.sz.sqrt_rc;             
%         
%     tmp0 = tmp0 + abs( V( :, :, pp )) .^ 2;
%     
% end
% 
% % need the sqrt of this when enforcing measurement constraint:
% tmp0 = sqrt( tmp0 );


V = fft( fftshift( fft( fftshift( phi, 1 ), [], 1 ), 2 ), [], 2 ) / sol.sz.sqrt_rc;
tmp0 = sqrt( sum( abs( V ) .^ 2, 3 ));

% % repmat2 = repmat( tmp0, [ 1, 1, sol.probe.scpm.N] );
% 
% norm( Btmp0 - tmp0, 'fro' )
% norm( V(:) - BV(:), 'fro' )

%{

for pp = 1 : sol.probe.scpm.N

    figure; imagesc( log10( 1 + abs( V( :, :, pp ))))


end

figure; imagesc( log10( 1 + tmp0 ))
figure; imagesc( log10( 1 + meas.D ))

%}

%==================================================================================================

% add option for fresnel propagator here

% exit waves corresponding to the different probe modes:
for pp = 1 : sol.probe.scpm.N             % size( phi, 3 )
  
    %==============
    
%     tmp1 = meas.D .* (( V( :, :, pp ) ./ ( 1e-7 + tmp0 )) .* meas.Dneq0 ) + V( :, :, pp ) .* meas.Deq0;
    tmp1 = meas .* ( V( :, :, pp ) ./ ( 1e-7 + tmp0 )) + V( :, :, pp ) .* meas_Deq0;
    
    %==============
    
    V( :, :, pp ) = fftshift( ifft2( measLPF .* tmp1 )) * sol.sz.sqrt_rc;
    
%     V( :, :, pp ) = fftshift( ifft2( measLPF .* meas.D .* ( V( :, :, pp ) ./ ( 1e-7 + tmp0 )))) * sol.sz.sqrt_rc;    

    %==============
    
end


%{

figure; 
imagesc(log10(1+abs(fftshift( meas.D .* ( V( :, :, pp ) ./ ( 1e-7 + tmp0 )) + V( :, :, pp ) .* meas.Deq0  )))); 
daspect([1 1 1]);
set(gca,'FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
colormap jet

%}

%==================================================================================================