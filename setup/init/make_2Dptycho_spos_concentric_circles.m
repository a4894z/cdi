function [ spos ] = make_2Dptycho_spos_concentric_circles( expt )

%================================================
% load a set of previously defined scan positions
%================================================

Z = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/rPIE_vs_MB_mat/no_noise/sim_ptycho2DTPA.mat', 'expt' );

spos = Z.expt.spos;

spos.indxsubset = spos.indx;

spos.rs0 = spos.rs;

% spos.rs = spos.rs0;

%=================================================================
% introduce misc probe position goofiness and scan position errors
%=================================================================

% translation
%
% spos.rs( :, 1 ) = spos.rs( :, 1 ) + 6;
% spos.rs( :, 2 ) = spos.rs( :, 2 ) - 20;

%=======================

% % shear
% 
% % sx = +0.5 * sign( 2 * rand - 1 );
% % sy = +0.5 * sign( 2 * rand - 1 );
% % sx = +0.5 * ( 2 * rand - 1 );
% % sy = +0.5 * ( 2 * rand - 1 );
% sx = +0.3;
% sy = -0.1;
% 
% spos.rs = transpose( [ 1, 0; sx, 1 ] * transpose( spos.rs ));
% 
% figure; 
% plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-400, 400])
% ylim([-400, 400])
% daspect([1 1 1])  
% grid on
%         
% spos.rs = transpose( [ 1, sy; 0, 1 ] * transpose( spos.rs ));
% 
% figure; 
% plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-400, 400])
% ylim([-400, 400])
% daspect([1 1 1])  
% grid on

%=======================

% rotate

figure; 
plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-600, 600])
ylim([-600, 600])
daspect([1 1 1])  
grid on



th = 12 * pi / 180;
% spos.rs = spos.rs * [ cos( th ), sin( th ); -sin( th ), cos( th ) ];
% spos.rs = fliplr( fliplr( spos.rs ) * [ cos( th ), +sin( th ); -sin( th ), cos( th ) ] );
spos.rs = spos.rs * [ cos( th ), -sin( th ); +sin( th ), cos( th ) ];

figure; 
plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-600, 600])
ylim([-600, 600])
daspect([1 1 1])  
grid on

%=======================

% scale
% 
% sx = 1.0 + 0.1;
% sy = 1.0 - 0.2;
% spos.rs = [ sy, 0; 0, sx ] * transpose( spos.rs );
% 
% 
% figure; 
% plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-600, 600])
% ylim([-600, 600])
% daspect([1 1 1])  
% grid on
% 
% 5;

%=======================
% random position errors
%=======================
% 
% spos.rs = spos.rs + round( 20 * ( 2 * rand( spos.N, 2 ) - 1 ));

%====================================================================================================================================================

return

%============================================================================
% create a 2 x ( # of scan pos ) array of scan positions with units of pixels
% define downstream optical axis as +y, then +x ( horizontal ) is inboard 
% (towards storage ring as viewed from beamline ), and +z (vertical) is up.
%============================================================================

% if isfield( expt, 'spos' ), spos = expt.spos; end

spos.type = 'concentric_circles';

spos.FOVx = 1.1 * 320;                       % fov in the x (col) direction in pixels
spos.FOVy = 1.7 * 320;                       % fov in the y (row) direction in pixels
spos.first_shell = 6;                  % number of scan locations in first circle
spos.shell_dr = 20;                    % shell width in pixels
spos.theta_offset = 0 * 20 * pi / 180;

%====================================================================================================================================================

spos.rs = [];

spos.rs( 1, : ) = [ 0, 0 ];

% Prepare field of view
rmax = sqrt( ( spos.FOVx / 2 ) ^ 2 + ( spos.FOVy / 2 ) ^ 2 );
nr = 1 + floor( rmax / spos.shell_dr );

for ir = 1 : ( nr + 1 )
  
  rr = ir * spos.shell_dr;
  dth = -2 * pi / ( spos.first_shell * ir );
  
  for ith = 0 : spos.first_shell * ir - 1
    
    th = ith * dth + spos.theta_offset;
    x2 = rr * cos( th );
    x1 = rr * sin( th );
    
    % if outside of desired fov, don't scan here
    if ( abs( x1 ) > spos.FOVy / 2 || ( abs( x2 ) > spos.FOVx / 2 ))
        continue; end
   
%     spos.rs( end + 1, : ) = round( [ x1, x2 ]); 
    spos.rs( end + 1, : ) = single( [ x1, x2 ] ); 
    
  end
  
end

spos.N = size( spos.rs, 1 );

spos.indx       = 1 : spos.N;
spos.indxsubset = spos.indx;






spos.rs0 = spos.rs;

% spos.rs = spos.rs0;


% figure; 
% plot_2Dscan_positions( spos.rs0, [], spos.rs( end, : ), [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-400, 400])
% ylim([-400, 400])
% daspect([1 1 1])  
% grid on



%=================================================================
% introduce misc probe position goofiness and scan position errors
%=================================================================

% translation
%
% spos.rs( :, 1 ) = spos.rs( :, 1 ) + 6;
% spos.rs( :, 2 ) = spos.rs( :, 2 ) - 20;

%=======================

% shear

% sx = +0.5 * sign( 2 * rand - 1 );
% sy = +0.5 * sign( 2 * rand - 1 );
% sx = +0.5 * ( 2 * rand - 1 );
% sy = +0.5 * ( 2 * rand - 1 );
sx = +0.1;
sy = -0.2;

spos.rs = spos.rs * [ 1, sx; 0, 1 ];

figure; 
plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-400, 400])
ylim([-400, 400])
daspect([1 1 1])  
grid on
        
spos.rs = spos.rs * [ 1, 0; sy, 1 ];

figure; 
plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-400, 400])
ylim([-400, 400])
daspect([1 1 1])  
grid on

%=======================

% % rotate
% 
% th = 30 * pi / 180;
% spos.rs = spos.rs * [ cos( th ), sin( th ); -sin( th ), cos( th ) ];

%=======================

% scale

sx = 1.0 + 0.1;
sy = 1.0 - 0.2;
spos.rs = spos.rs * [ sy, 0; 0, sx ];


figure; 
plot_2Dscan_positions( spos.rs0, [], spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-600, 600])
ylim([-600, 600])
daspect([1 1 1])  
grid on

5;

%=======================

% randomize positions

% spos.rs = spos.rs + round(  20 * ( 2 * rand( spos.N, 2 ) - 1 ));

% spos.rs( :, 1 ) = spos.rs( :, 1 ) + round( spos.randR * ( 2 * rand( spos.N, 1 ) - 1 ));
% spos.rs( :, 2 ) = spos.rs( :, 2 ) + round( spos.randC * ( 2 * rand( spos.N, 1 ) - 1 ));

%====================================================================================================================================================





    
    
    
    
    
    
    
    
    
    
    