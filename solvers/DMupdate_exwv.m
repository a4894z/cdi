function [ phi ] = DMupdate_exwv( phi, sol, expt )
% DMmeasupdate2DTPA_allspos

% j is iteration index
% s is index for scan pos
% p is index for probe modes
% t is index for sample transfer function modes

% phi_{ j+1, s, p } <-- phi_{ j, s, p } + beta ( proj_M[ 2 P_{ j, p } T_{ j, s } - phi_{ j, s, p } ] - P_{ j, p } T_{ j, s } )

%==================================================================================================

% Ps_phi = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );    
% Rs_phi = zeros( sol.sz.r, sol.sz.c, sol.probe.scpm.N, 'single' );    

%==================================================================================================

for ss = sol.spos.updateorder
  
    rs = sol.spos.rs( ss, : );          % !!!!!!!!!!!
    
    [ Ps_phi, ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, rs, sol.spos.shifttype );
    
    Rs_phi = 2 * Ps_phi - phi( :, :, :, ss );

    %==============
    
    % enforce measurement on exitwave view reflection operator:
    % TODO: propagator_type = { 'fft', 'frsnl', 'frsnl_23', 'frsnl_213' }, USING FUNCTION HANDLE
    Pm_Rs_phi = enforce_2DTPAmeas( Rs_phi, expt.meas.SI( ss ), sol.measLPF, sol );

    %==============
    
    % difference map update:

    phi( :, :, :, ss ) = phi( :, :, :, ss ) + sol.betaDM * ( Pm_Rs_phi - Ps_phi ); 

    %==============

end



