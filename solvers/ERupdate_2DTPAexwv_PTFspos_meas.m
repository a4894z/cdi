function [ phi ] = ERupdate_2DTPAexwv_PTFspos_meas( sol, expt )

% projapprx <--> PA
% multslice <--> MS

% ePIEexwvupdt_3DBprojapprx    % 3D sample, bragg, projection approx
% ePIEexwvupdt_3DTmultslice    % 3D sample, transmission, multislice
% ePIEexwvupdt_2DRprojapprx    % 2D sample, reflection, projection approx
% ePIEexwvupdt_2DTprojapprx    % 2D sample, transmission, projection approx 

%==================================================================================================

% j is iteration index
% s is index for scan pos
% p is index for probe modes
% t is index for sample transfer function modes

% difference map:
% phi_{j+1, s, p} <-- phi_{j, s, p} + beta ( proj_M[ 2 P_{j, p} T_{j, s} - phi_{j, s, p} ] - P_{j, p} T_{j, s} )

% measurement projection operator:
% proj_M[ phi_{s,p} ] = F^{-1}[ measSI_{s} * F[ phi_{s,p} ] / sqrt( sum_p | F[ phi_{s,p} ] |^2 ) ]

%==================================================================================================

phi = zeros( sol.sz.sz( 1 ), sol.sz.sz( 2 ), sol.probe.scpm.N, sol.spos.N, 'single' );      % probe .* sample view 

%==================================================================================================

% if ~isfield( sol.spos, 'updateorder' )
%     
%     updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));
%     
% else
%     
%     updateorder = sol.spos.updateorder;
%     
% end

%==================================================================================================

for ss = 1 : sol.spos.N   % order * DOESN'T * matter here !!!
% for ss = 1 : length( sol.spos.indxsubset )   % order * DOESN'T * matter here !!!
% for ss = updateorder  % order * DOESN'T * matter here !!!
    
%     rs = +sol.spos.rs( ss, : );    % !!!!!!!!!!!
%     rs = -sol.spos.rs( ss, : );
    
%     rs = [ +rs(1), -rs(2) ];
%     rs = [ -rs(1), +rs(2) ];
%     rs = [ -rs(1), -rs(2) ];
    
%     rs = -sol.spos.rs( ss, : );
    
    %===========================
    
    % form 2D exitwave(s) under projection approximation in transmission geometry:
    [ phi( :, :, :, ss ), ~ ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, sol.spos.rs( ss, : ), sol.spos.shifttype );
    
    %===========================
  
    % opts.z2z3prop = { 'fftn', 'fft2', 'frsnl_23', 'frsnl_213' }
    
    % enforce measurement on exit waves:
    phi( :, :, :, ss ) = enforce_2DTPAmeas( phi( :, :, :, ss ), expt.meas.SI( ss ), sol.measLPF, sol );

    %===========================

end
