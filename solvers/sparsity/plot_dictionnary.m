function H = plot_dictionnary( DL, nb )

D = DL.D;

I = round( linspace( 1, size( D, 2 ), prod( nb )));

pad = 5;

H = complex( zeros( nb( 1 ) * ( DL.wy + pad ), nb( 2 ) * ( DL.wx + pad ), 'single'));

k = 0;

for i = 1 : nb( 2 )
    for j = 1 : nb( 1 )
        
        k = k + 1;
        
        if k <= length( I ) 
            
            v = D( :, I( k ));
            
            v = reshape( v, DL.wy, DL.wx, 1 );

            selx = ( ( i - 1 ) * ( DL.wx + pad ) + 1 ) : ( ( i - 1 ) * ( DL.wx + pad ) + DL.wx );
            sely = ( ( j - 1 ) * ( DL.wy + pad ) + 1 ) : ( ( j - 1 ) * ( DL.wy + pad ) + DL.wy );

            H( sely, selx, : ) = v;
                
        end
    end
end

H( end - pad + 1 : end, :, : ) = [];
H( :, end - pad + 1 : end, : ) = [];
% H( H == 0 ) = 10;


subplot(131); 
imagescHSV( log10( 1 + 10^2 * abs( H )) .* exp( 1i * angle( H ))); 
% imagescHSV( H ); 
title('hsv'); 
daspect([1 1 1]);

ax1 = subplot(132); 
imagesc( real( H )); 
title('re'); 
daspect([1 1 1]); 
colormap( ax1, bone )

ax2 = subplot(133); 
imagesc( imag( H )); 
title('im'); 
daspect([1 1 1])
colormap( ax2, bone )




% subplot(131); 
% imagescHSV( log10( 1 + 10^-0 * abs( H )) .* exp( 1i * angle( H ))); 
% % imagescHSV( H ); 
% title('hsv'); 
% daspect([1 1 1]);
% 
% ax1 = subplot(132); 
% imagesc( log10( 1 + 10^-0 * abs( H ))); 
% % imagesc( abs( H )); 
% title('abs'); 
% daspect([1 1 1]); 
% colormap( ax1, bone )
% 
% ax2 = subplot(133); 
% imagesc( angle( H )); 
% title('angle'); 
% daspect([1 1 1])
% colormap( ax2, hsv )




% data( ~any(data,2), : ) = [];  % remove zero rows
% data( :, ~any(data,1) ) = [];  % remove zero columns






%{

D = DL.D;

I = round( linspace( 1, size( D, 2 ), prod( nb )));

% a = 1 : size( D, 2 );
% numelements = prod( nb );
% indices = randperm( length( a ));
% indices = indices( 1 : numelements );
% I = a(indices);

pad = 5;

H = complex( zeros( nb( 1 ) * ( DL.wy + pad ), nb( 2 ) * ( DL.wx + pad ), 'single'));

k = 0;

for i = 1 : nb( 2 )
    for j = 1 : nb( 1 )
        
        k = k + 1;
        
        if k <= length( I ) 
            
            v = D( :, I( k ));
            
            v = reshape( v, DL.wy, DL.wx, 1 );

            selx = ( ( i - 1 ) * ( DL.wx + pad ) + 1 ) : ( ( i - 1 ) * ( DL.wx + pad ) + DL.wx );
            sely = ( ( j - 1 ) * ( DL.wy + pad ) + 1 ) : ( ( j - 1 ) * ( DL.wy + pad ) + DL.wy );

            H( sely, selx, : ) = v;
                
        end
    end
end

H( end - pad + 1 : end, :, : ) = [];
H( :, end - pad + 1 : end, : ) = [];
% H( H == 0 ) = 10;

subplot(131); 
imagescHSV( log10( 1 + 10^-0 * abs( H )) .* exp( 1i * angle( H ))); 
% imagescHSV( H ); 
title('hsv'); 
daspect([1 1 1]);

ax1 = subplot(132); 
imagesc( log10( 1 + 10^-0 * abs( H ))); 
% imagesc( abs( H )); 
title('abs'); 
daspect([1 1 1]); 
colormap( ax1, gray )

ax2 = subplot(133); 
imagesc( angle( H )); 
title('angle'); 
daspect([1 1 1])
colormap( ax2, hsv )



    
%}
