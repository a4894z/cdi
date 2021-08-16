function [ Cs, lam, phi ] = ADMMupdate_exwvmeas_sDLsparsecode( Cs, lam, meas, sol )

Vs = Cs - lam;



abs_C = abs( Vs );
thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * sDL.sparse_lvl );

Vs = Vs .* ( abs_C >= thresh );                                                             % hard thresh
% Vs = ( abs_C -  1.0 * thresh ) .* ( abs_C >= thresh ) .* exp( 1i * angle( Vs ));          % soft thresh




Vm = Vs + lam;


[ phi ] = trainingsetpatches2image( sDL.D * Vm, sDL, sol.sz.sz );

phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );


[ Ip, ~ ] = image2trainingsetpatches( phi, sDL, sol.sz.sz );

Cs = sDL.D' * Ip;







lam = lam + Vs - Cs;