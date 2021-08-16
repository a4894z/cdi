function [ Vm, Vs, lam, sDL ] = ADMMupdate_exwv_DsCs_meas( Vm, lam, sDL, meas, sol )

    %==============================================================================================
    
    Vs = Vm - lam;
    
    [ Ip_phi, ~ ] = image2trainingsetpatches( Vs, sDL, sol.sz.sz );

    sDL.D = sDL_dictionary_update( Ip_phi, sDL );
    sDL.C = sDL_sparsecode_update( Ip_phi, sDL );

    
%     Cs = sDL_sparsecode_update( Ip_phi, sDL );
%     Ds = sDL_dictionary_update( Ip_phi, sDL );
%     
%     sDL.C = Cs;
%     sDL.D = Ds;
    
    [ Vs ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sol.sz.sz );
    
    %==============================================================================================

    Vm = Vs + lam;

    Vm = enforce_2DTPAmeas( Vm, meas, sol.measLPF, sol );
    
    %==============================================================================================

    lam = lam + 0.9 * ( Vs - Vm );
