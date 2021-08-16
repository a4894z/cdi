function [ V ] = make_2Dgaussian( N, mu, sdev )

    %==================

%     N    = single( N ); 
%     mu   = single( mu );
%     sdev = single( sdev );
    
    [ c, r ] = meshgrid( 1 : N( 2 ), 1 : N( 1 ) );

%     c = single( c ); 
%     r = single( r ); 

    %==================

    V_c = exp( -1 * ( c - mu( 2 ) ).^2 / ( 1e-7 + 2 * sdev( 2 ) )^2 );
    V_r = exp( -1 * ( r - mu( 1 ) ).^2 / ( 1e-7 + 2 * sdev( 1 ) )^2 );

    V = V_r .* V_c;

    %==================