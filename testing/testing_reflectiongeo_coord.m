%





clear; close all

Nr = 512; 
Nc = 512;
z23 = 3.0;
kev = 7.0;
llambda = ( 12.4 / kev ) * 1e-10;

% pixel size at measurement plane z3:
csys.z3.dLy = 55e-6;                            % ( in meters )
csys.z3.dLx = 55e-6;                            % ( in meters ) 

cols = (( -0.5 * Nc + 1 ) : ( 0.5 * Nc )) * csys.z3.dLx;
rows = (( -0.5 * Nr + 1 ) : ( 0.5 * Nr )) * csys.z3.dLy;
% cols = (1 : Nc ) * csys.z3.dLx;
% rows = (1 : Nr ) * csys.z3.dLy;


[ y, z ] = meshgrid( cols, rows );




alpha_i = 0.5 * pi / 180;
alpha_f   = asin( z ./ sqrt(z23 ^ 2 + y .^ 2 + z .^ 2 ));
two_theta = atan(  y / z23 );



figure;
imagesc( cols / csys.z3.dLx, rows / csys.z3.dLy, two_theta )
colormap jet
colorbar
xlabel('column pixel index')
ylabel('row pixel index')
title('2 theta = atan( y / d_1 )')




figure;
imagesc( cols / csys.z3.dLx, rows / csys.z3.dLy, alpha_f  )
colormap jet
colorbar
xlabel('column pixel index')
ylabel('row pixel index')
title('alpha_f = asin( z / sqrt{d_1^2 + y^2 + z^2 })')





q_x = 1e-9 * ( 2 * pi / llambda ) * ( cos( alpha_f ) .* cos( two_theta ) - cos( alpha_i ));

figure; imagesc( q_x  ); colormap jet; colorbar
title('q_x')



q_y = 1e-9 * ( 2 * pi / llambda ) * ( cos( alpha_f ) .* sin( two_theta )) ;

figure; imagesc( q_y ); colormap jet; colorbar
title('q_y')



q_z = 1e-9 * ( 2 * pi / llambda ) * ( sin( alpha_f ) + sin( alpha_i ));

figure; imagesc( q_z ); colormap jet; colorbar
title('q_z')











 [ V3 ] = make_2Dellipsoid( [Nr, Nc], 0.3 * [Nr, Nc] );









