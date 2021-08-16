
%{

clear; close all;
% restoredefaultpath; 
expt.paths.code = pwd;
addpath( genpath( expt.paths.code ));
clearvars -except expt sol

rng( 'shuffle' )
try gpuDevice, reset( gpuDevice( 1 )); catch, end;


clear; close all; testing_ndim_indexing




% clear; close all
% 
% indx = [[ 7; 13 ], [ 22; 32 ], [ 51; 60 ]];
% a    = [[ 1; 2 ], [ 3; 4 ], [ 5; 6] ]
% 
% b = zeros( 5, 4, 3  );
% 
% b( indx ) = a

% A = 1;
% B = 1;
% C = 1; 
% D = 1;
% [x y] = meshgrid(-1:0.1:1); % Generate x and y data
% z = -1/C*(A*x + B*y + D); % Solve for z data
% surf(x,y,z) %Plot the surface

%}
%==================================================================================================

clear; close all; 
paths_code = '~/Documents/Science/Matlab/Code/tomolamino/testing';
cd( paths_code )

restoredefaultpath; 
addpath( genpath( '~/Documents/Science/Matlab/Code/cdi' ));
addpath( genpath( paths_code ));
clearvars -except expt sol paths_code

%==================================================================================================
% 2d

clear; close all;

Nr = 64;
Nc = 64;

lnwidth = 2.0;

[ cc, rr ] = meshgrid( 1 : Nc, 1 : Nr );




[ X, Y, Z ] = meshgrid( -5 : 0.2 : 5 );
V = X .* exp( -X .^ 2 - Y .^ 2 - Z .^ 2 );

[ xsurf, ysurf ] = meshgrid( -2 : 0.2 : 2 );
zsurf = xsurf .^ 2 - ysurf .^ 2;
slice( X, Y, Z, V, xsurf, ysurf, zsurf )








r1a = ((5.5 * cc + 0 * lnwidth - 32 ) < rr);
r1b = ((5.5 * cc + 1 * lnwidth - 32 ) < rr);
r1 = ( r1a - r1b );

figure; 
imagesc( 1.0 * r1 )










r1a = ((5.5 * cc + 0 * lnwidth - 32 ) < rr);
r1b = ((5.5 * cc + 1 * lnwidth - 32 ) < rr);
r1 = ( r1a - r1b );

r2a = ((2.5 * cc + 1 * lnwidth ) < rr);
r2b = ((2.5 * cc + 2 * lnwidth ) < rr);
r2 = ( r2a - r2b );

r3a = ((2.5 * cc + 2 * lnwidth ) < rr);
r3b = ((2.5 * cc + 3 * lnwidth ) < rr);
r3 = ( r3a - r3b );

figure; 
imagesc( 1.0 * r1 + 2.0 * r2 + 3.0 * r3)



return

%==================================================================================================
% 3d





[ cc, rr, pp ] = meshgrid( single( 1 : 256 ), single( 1 : 256 ), single( 1 : 512 ));

Vbp1 = ( rr < (( 0.5 * cc + 0.4 * pp ) - 0 ));

Vbp2 = ( rr <= (( 0.5 * cc + 0.4 * pp ) - 0 ));

figure; isosurface( Vbp1 )
figure; isosurface( Vbp2 )
figure; isosurface( Vbp1 - Vbp2 )



V = Vbp1 - Vbp2;
for pp = 1 : 512
    
    figure( 555 ); imagesc( V( :, :, pp ))

 
end


Vbp1 = ((rr - (( 0.5 * cc + 0.4 * pp ) - 0 )) < 1.0);


for pp = 1 : 512
    
    figure( 555 ); imagesc( Vbp1( :, :, pp ))

 
end



return

%==================================================================================================

Nr = 5;
Nc = 4;

A0v = 10 * transpose( 1 : ( Nr * Nc )) + 0

A0m = reshape( A0v, [ Nr, Nc ] );
Am  = zeros( Nr, Nc,  3 );

ind( :, 1 ) = [ 24, 25, 26, ...
                34, 35, 36, ...
                44, 45, 46, ...
                54, 55, 56 ];

ind( :, 2 ) = [ 54, 55, 56, ...
                64, 75, 76, ...
                74, 75, 76, ...
                04, 05, 06 ];


ind( :, 3 ) = [ 42, 43, 44, ...
                52, 53, 54, ...
                62, 63, 64, ...
                72, 73, 74 ];


A0m( ind );
A0v( ind )

tmp0 = A0m( ind );
% tmp0 = reshape( tmp0, [ 3, 4, 3] )


% Av( ind, : ) = tmp0

Av  = zeros( Nr * Nc, 3 );

% Av( ind ) = tmp0

Av( ind( 1, :) ) = tmp0( 1, : )







Av  = zeros( Nr * Nc * 1, 1 );
Av( ind ) = tmp0


reshape( Av, [ Nr * Nc, 3 ] )







Av( ind( :, 1:2 ), 1:2) = tmp0( :, 1:2 )




Av( ind( :, 1 ), 1 ) = tmp0( :, 1 )
Av( ind( :, 2 ), 2 ) = tmp0( :, 2 )
Av( ind( :, 3 ), 3 ) = tmp0( :, 3 )






Am( :, :, 1 ) = reshape( Av(:, 1), [Nr, Nc] )

%==================================================================================================


clear; close all;

B = [ 1, 2, 3, 4; 
      9, 8, 7, 6; 
      4, 5, 6, 7;
      6, 5, 4, 3;
      7, 8, 9, 10;
      3, 4, 5, 6;  ];

A = transpose( 1 : 10 ); 
A = repmat( A, 1, 5 );

for cc = 2 : size( A, 2 )
    
    A( :, cc ) = 10 * A( :, cc - 1 );

end

% A( :, :, 2 ) = A( :, : );
% A( :, 1, 2 ) = 10 * A( :, 1, 2 ) + 1;
% 
% 
% for cc = 2 : size( A, 2 )
%     
%     A( :, cc, 2 ) = 10 * A( :, cc - 1, 2 );
% 
% end











C = A( B, : )




E = reshape( C, [ size( B ), size( C, 2 ) ] )



D = A( B, :, : );

F = reshape( D, [ size( B ), size( D, 2 ), size( D, 3 ) ] )










clear; close all;


z = [ [ 11; 21; 31 ], [ 101; 201; 301 ], [ 701; 801; 901 ], 60 * [ rand; rand; rand ] ]; 
v = [ [  1;  2;  3 ], [   1;   3;   4 ], [   4;   1;   3 ],      [    2;    1;    4 ] ] + [0, 4, 8, 12]; 

b = [ [ 0, 0, 0, 0 ]; [ 0, 0, 0, 0 ]; [ 0, 0, 0, 0 ]; [ 0, 0, 0, 0 ] ]'; 




b( v ) = z





















