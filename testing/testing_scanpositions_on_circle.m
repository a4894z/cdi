


close all

y1 = 4 * ( 2*rand - 1 ); 
x1 = 3 * ( 2*rand - 1 ); 
y1c = 24 *( 2*rand - 1 ); 
x1c = 33 *( 2*rand - 1 ); 

R1 = sqrt( y1^2 + x1^2 ); 
R2 = ( 0.1 * rand + 0.7 ) * R1; 

th = 0 : pi/1000 : 2*pi;

x2 = linspace( -7, 7, length( th ));

sinth = y1 / R1;

xunit = R1 * cos(th) + x1c;
yunit = R1 * sin(th) + y1c;

figure;
grid on

hold on
h = plot(xunit, yunit);
hold off

xunit = R2 * cos(th) + x1c;
yunit = R2 * sin(th) + y1c;

hold on
h = plot(xunit, yunit);
hold off

hold on
plot(x1,y1,'r*')
hold off

hold on
plot(x2,y1/x1 * x2)
hold off
daspect([1 1 1])

y2 = y1 * R2 / R1
x2 = x1 * R2 / R1






tmp0 = expt.spos.rs( expt.spos.indxsubset, : );

R1 = sqrt(( spos_rs( :, 1 ) - tmp0( :, 1 )).^2 + ( spos_rs( :, 2 ) - tmp0( :, 2 )).^2 );
R1 = sqrt(( spos_rs( :, 1 ) - tmp0( :, 1 )).^2 + ( spos_rs( :, 2 ) - tmp0( :, 2 )).^2 );
R2 = 12;
Ireset = ( R2 < R1);

spos_rs( Ireset, 1 ) = tmp0( Ireset, : );

y2 = y1 * R2 / R1
x2 = x1 * R2 / R1
