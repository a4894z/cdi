function [ phi, sDL ] = ERupdate_exwv_sDL_C_meas( phi, sDL, meas, sol )  % ( phi, sDL, meas, sol, niter )


%==================================================================================================
% EXIT WAVE MEASUREMENT CONSTRAINT UPDATE:

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );

% tmp0 = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
% phi = 2 * tmp0 - phi;

% phi = phi .* sol.probe.S;

%==================================================================================================
% SYNTHESIS SPARSE CODE UPDATE

[ Ip, sDL ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );

sDL.C = sDL_sparsecode_update( Ip, sDL );
% sDL.D = sDL_dictionary_update( Ip, sDL.C, sDL.D );

[ phi ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sol.sz.sz );

%==================================================================================================


