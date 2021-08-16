function [ D ] = dct2_basis_functions( wx, wy )



[ nn, mm ] = meshgrid( 0 : wx- 1, 0 :  wy - 1 );

nr_max = wy; 
nc_max = wx;

D = zeros( wy * wx, nr_max * nc_max, 'single' );
Ds = zeros( wy * wx, nr_max * nc_max, 'single' );

ii = 1;
jj = 1;

% for nr = 0 : nr_max - 1
%     for nc = 0 : nc_max - 1
        
for nr = 0 : nr_max - 1
    for nc = 0 : nc_max - 1

        
        if nr == 0, alpha_nr = 1 / sqrt( wy );
        else, alpha_nr = sqrt( 2 /  wy ); end

        if nc == 0, alpha_nc = 1 / sqrt( wx );
        else, alpha_nc = sqrt( 2 /  wx ); end
        
        B_pq = ( 0 + 1 * alpha_nr * alpha_nc ) * cos( pi * ( 2 * mm + 1 ) * nr / ( 2 * wy )) .* cos( pi * ( 2 * nn + 1 ) * nc / ( 2 * wx ));

%         B_pq_sin = ( 0 + 1 * alpha_nr * alpha_nc ) * sin( pi * ( 2 * mm + 1 ) * nr / ( 2 * wy )) .* sin( pi * ( 2 * nn + 1 ) * nc / ( 2 * wx ));
%         B_pq = B_pq + 1i * B_pq_sin;

%         B_pq = ( 0 + 1 * alpha_nr * alpha_nc ) * exp( 1i * pi * ( 2 * mm + 1 ) * nr / ( 2 * wy )) .* exp( 1i * pi * ( 2 * nn + 1 ) * nc / ( 2 * wx ));
        
        % B_pq = B_pq / sqrt( wy * wy );
        

%         B_pq = B_pq .* exp( 1i * B_pq );
        

        
        
        D( :, ( ii - 1 ) * nc_max + jj ) = B_pq( : );

%         Ds( :, ( ii - 1 ) * nc_max + jj ) = Bs_pq( : );
      
        jj = jj + 1;
        
    end
    
    ii = ii + 1;
    jj = 1;
    
end



% figure; imagesc(  abs( D' * D ))
% 5;
% 
% figure; imagesc(  abs( Ds' * Ds ))
% 5;