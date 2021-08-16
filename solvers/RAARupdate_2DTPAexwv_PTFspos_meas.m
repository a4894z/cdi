function [ phi ] = RAARupdate_2DTPAexwv_PTFspos_meas( phi, sol, expt )

% sol.spos.updateorder = randperm( length( sol.spos.indxsubset ));
% for ss = sol.spos.updateorder     

%==================================================================================================

if ~isfield( sol.spos, 'updateorder' )
    
    updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));
    
else
    
    updateorder = sol.spos.updateorder;
    
end

%==================================================================================================

for ss = updateorder   % order * DOESN'T * matter here !!!

    %==============

    [ Ps_phi, ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, +sol.spos.rs( ss, : ), sol.spos.shifttype );

    Rs_phi = 2 * Ps_phi - phi( :, :, :, ss );

    %==============
    
    % TODO: propagator_type = { 'fft', 'frsnl', 'fresnel23', 'fresnel213' }
    
    % enforce measurement on exitwave views:
%     Pm_phi = enforce_2DTPAmeas( phi( :, :, :, ss ), expt.meas.SI( ss ), sol.measLPF, sol );
    Pm_Ps_phi = enforce_2DTPAmeas( Ps_phi, expt.meas.D( :, :, ss ), expt.meas.Deq0( :, :, ss ), sol.measLPF, sol );
    
    % enforce measurement on exitwave view reflection operator:
    Pm_Rs_phi = enforce_2DTPAmeas( Rs_phi, expt.meas.SI( ss ), sol.measLPF, sol );
   
    Rm_Rs_phi = 2 * Pm_Rs_phi - Rs_phi;
    
    %==============
    
    % RAAR update:
%     phi( :, :, :, ss ) = 0.5 * sol.RAAR.beta * ( Rm_Rs_phi + phi( :, :, :, ss )) + ( 1 - sol.RAAR.beta ) * Pm_phi;
    phi( :, :, :, ss ) = 0.5 * sol.RAAR.beta * ( Rm_Rs_phi + phi( :, :, :, ss )) + ( 1 - sol.RAAR.beta ) * Pm_Ps_phi;
        
end


