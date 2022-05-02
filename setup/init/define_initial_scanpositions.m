function [ sol, expt ] = define_initial_scanpositions( sol, expt )

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

%==========================================
% Other transformation(s) on scan positions
%==========================================
% 
% % scaling in rows/cols
% sol.spos.rs( :, 1 ) = 0.95 * sol.spos.rs( :, 1 );
% sol.spos.rs( :, 2 ) = 1.0 * sol.spos.rs( :, 2 );

%==============================================================================================
% Remove some scan positions from the reconstruction, remove corresponding measurements as well
%==============================================================================================

% rs_skip = 2;
% 
% sol.spos.indxsubset = sol.spos.indxsubset( 1 : rs_skip : end );  
% sol.spos.rs         = sol.spos.rs( sol.spos.indxsubset, : );

%========
% 
% Irt = sol.spos.rs( :, 1 ) < +5000;
% Irb = sol.spos.rs( :, 1 ) > -5000;
% 
% Icl = sol.spos.rs( :, 2 ) < +2400;
% Icr = sol.spos.rs( :, 2 ) > -2400;
% 
% sol.spos.indxsubset = sol.spos.indxsubset( Icl & Icr & Irt & Irb );
% sol.spos.rs         = sol.spos.rs( sol.spos.indxsubset, : );
% 
% %================================
% % Recenter the modified positions
% %================================
% 
% sol.spos.rs = sol.spos.rs - min( sol.spos.rs, [], 1 );
% sol.spos.rs = sol.spos.rs - max( sol.spos.rs, [], 1 ) * 0.5;  
% 
% %============================
% % Get new number of positions
% %============================
% 
% sol.spos.N  = size( sol.spos.rs, 1 );
% 
% %================================================================
% % If we removed scan positions, remove corresponding measurements
% %================================================================
% 
% expt.meas.D    = expt.meas.D( :, :, sol.spos.indxsubset );
% expt.meas.Deq0 = ( expt.meas.D == 0 );

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

