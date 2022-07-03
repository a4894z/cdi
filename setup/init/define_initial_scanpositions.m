function [ sol, expt ] = define_initial_scanpositions( sol, expt )

%=======================
% FOR OLD SIMULATED DATA
%=======================

% try expt.spos.indxsubset = expt.spos.indx; catch, end
   
%=========================
% Load old scan positions?
%=========================

% old_spos = load( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/data/zjiang202112/probes/L0105_to_L0113_combined_512x512_12Jan2022_t130245_it15001.mat', 'spos' );
% sol.spos = old_spos.spos;

%===================================================
% Use scan positions from what the beamline tells us
%===================================================
 
if ~isfield( sol, 'spos' ) || sol.init.reinit_spos

    sol.spos.N          = expt.spos.N;
    sol.spos.rs         = expt.spos.rs;
    sol.spos.indxsubset = expt.spos.indxsubset;
    
end

%=============================================
% Shearing transformation(s) on scan positions
%=============================================

% s_x = -0.01;
% sol.spos.shear_x = [ [ 1,   0 ]; ...
%                      [ s_x, 1 ] ];
% 
% sol.spos.rs = transpose( sol.spos.shear_x * transpose( sol.spos.rs ));

% s_y = -0.01;
% sol.spos.shear_y = [ [ 1, s_y ]; ...
%                      [ 0, 1   ] ];
%                  
% sol.spos.rs = transpose( sol.spos.shear_y * transpose( sol.spos.rs ));

%=============================================
% Rotation transformation(s) on scan positions
%=============================================

% theta = -1.0;
% sol.spos.rotation = [ [ +cosd( theta ), +sind( theta )  ]; ...
%                       [ -sind( theta ), +cosd( theta )  ] ];
% 
% sol.spos.rs = transpose( sol.spos.rotation * transpose( sol.spos.rs ));

% figure; 
% plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% daspect([1 1 1])  

%==========================================
% Other transformation(s) on scan positions
%==========================================

% % scaling in rows/cols
% sol.spos.rs( :, 1 ) = 1.00 * sol.spos.rs( :, 1 );
% sol.spos.rs( :, 2 ) = 1.00 * sol.spos.rs( :, 2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Remove some scan positions from the reconstruction, remove corresponding measurements as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================
% skip every other using a stride
%================================

% rs_skip = 2;
% 
% % sol.spos.indxsubset = sol.spos.indxsubset( 2 : rs_skip : end );         % even positions
% sol.spos.indxsubset = sol.spos.indxsubset( 1 : rs_skip : end );         % odd positions
% sol.spos.rs         = sol.spos.rs( sol.spos.indxsubset, : );

%=========================
% randomly choose a subset
%=========================

% sol.spos.indxsubset = single( randperm( sol.spos.N, round( 0.50 * sol.spos.N ) ));
% sol.spos.rs         = sol.spos.rs( sol.spos.indxsubset, : );

%==============================
% choose within a parallelogram
%==============================

%{

figure; 
plot_2Dscan_positions( sol.spos.rs, [], sol.spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  



            theta = +35.0;
            spos_rotation = [ [ +cosd( theta ), +sind( theta )  ]; ...
                              [ -sind( theta ), +cosd( theta )  ] ];
            
            rs_rot = transpose( spos_rotation * transpose( sol.spos.rs ));

figure; 
plot_2Dscan_positions( sol.spos.rs, [], rs_rot, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  

            
            

            Irt = ( rs_rot( :, 1 ) < +140 );
            Irb = ( rs_rot( :, 1 ) > -180 );
            
            Icl = ( rs_rot( :, 2 ) < +650 );
            Icr = ( rs_rot( :, 2 ) > -650 );
            
            sol.spos.indxsubset = sol.spos.indxsubset( Icl & Icr & Irt & Irb );
            rs_rot_subset1      = rs_rot( sol.spos.indxsubset, : );


figure; 
plot_2Dscan_positions( sol.spos.rs, [], rs_rot_subset1, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  
  

            theta = -1 * theta;
            spos_rotation = [ [ +cosd( theta ), +sind( theta )  ]; ...
                              [ -sind( theta ), +cosd( theta )  ] ];
               
            rs_rot_subset2 = transpose( spos_rotation * transpose( rs_rot_subset1 ));


figure; 
plot_2Dscan_positions( sol.spos.rs, [], rs_rot_subset2, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  



sol.spos.rs = rs_rot_subset2;
            
%}

%=======================================
% choose within a box using inequalities
%=======================================

% Irt = ( sol.spos.rs( :, 1 ) < +600 );
% Irb = ( sol.spos.rs( :, 1 ) > -400 );
% 
% Icl = ( sol.spos.rs( :, 2 ) < +8000 );
% Icr = ( sol.spos.rs( :, 2 ) > -8000 );
% 
% sol.spos.indxsubset = sol.spos.indxsubset( Icl & Icr & Irt & Irb );
% sol.spos.rs         = sol.spos.rs( sol.spos.indxsubset, : );

%================================
% Recenter the modified positions
%================================

% sol.spos.rs = sol.spos.rs - min( sol.spos.rs, [], 1 );
% sol.spos.rs = sol.spos.rs - max( sol.spos.rs, [], 1 ) * 0.5;  

%============================
% Get new number of positions
%============================

sol.spos.N = size( sol.spos.rs, 1 );

%========================================================================================
% If we removed scan positions, update corresponding used measurements and scan positions
%========================================================================================

expt.meas.D    = expt.meas.D( :, :, sol.spos.indxsubset );
expt.meas.Deq0 = ( expt.meas.D == 0 );

%========

sol.spos.rs0 = expt.spos.rs( sol.spos.indxsubset, : );





%{

figure; 
plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
daspect([1 1 1])  

% tmp0 = (max(sol.spos.rs(:,2)) - min(sol.spos.rs(:,2)));
% 
% hold on
% plot(( 0 + 0.1 * linspace( min(sol.spos.rs(:,2)), max(sol.spos.rs(:,2)), 1.0*(max(sol.spos.rs(:,2)) - min(sol.spos.rs(:,2))) )))
% hold off

%}

