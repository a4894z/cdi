
%
%{

clear; close all; testing_mb_vs_rPIEalpha

%}

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/';

%====================================================================================================================================================

%===============
% Stoch Block GD
%===============

path_data_stoch_bgd = {};
N_trials_stoch_bgd  = [];

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha1p00_oldsample_for_probe_update/' ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p10_oldsample_for_probe_update/' ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, 'blind_block_stoch/cdi_rPIE_stochGD_alpha0p01_oldsample_for_probe_update/' ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

path_data_stoch_bgd{ end + 1 } = [ rootpath_data, ];
N_trials_stoch_bgd( end + 1 )  = 10;

%=============
% 1% Minibatch
%=============

path_data_0p01 = {};
N_trials_0p01  = [];

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, 'mb0p01/cdi_rPIE_mb0p01_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, ];
N_trials_0p01( end + 1 )  = 10;

path_data_0p01{ end + 1 } = [ rootpath_data, ];
N_trials_0p01( end + 1 )  = 10;

%=============
% 5% Minibatch
%=============

path_data_0p05 = {};
N_trials_0p05  = [];

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p10_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p00000001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

path_data_0p05{ end + 1 } = [ rootpath_data, 'mb0p05/cdi_rPIE_mb0p05_alpha0p000000001_oldsample_for_probe_update/' ];
N_trials_0p05( end + 1 )  = 10;

%==============
% 10% Minibatch
%==============

path_data_0p10 = {};
N_trials_0p10  = [];

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p00000001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

path_data_0p10{ end + 1 } = [ rootpath_data, 'mb0p10/cdi_rPIE_mb0p10_alpha0p000000001_oldsample_for_probe_update/' ];
N_trials_0p10( end + 1 )  = 10;

%==============
% 20% Minibatch
%==============

path_data_0p20 = {};
N_trials_0p20  = [];

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p10_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p00000001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

path_data_0p20{ end + 1 } = [ rootpath_data, 'mb0p20/cdi_rPIE_mb0p20_alpha0p000000001_oldsample_for_probe_update/' ];
N_trials_0p20( end + 1 )  = 10;

%==============
% 33% Minibatch
%==============

path_data_0p33 = {};
N_trials_0p33  = [];

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p10_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p00000001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

path_data_0p33{ end + 1 } = [ rootpath_data, 'mb0p33/cdi_rPIE_mb0p33_alpha0p000000001_oldsample_for_probe_update/' ];
N_trials_0p33( end + 1 )  = 10;

%==============
% 50% Minibatch
%==============

path_data_0p50 = {};
N_trials_0p50  = [];

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha1p00_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p10_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p01_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p0001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p00001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p000001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p0000001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p00000001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

path_data_0p50{ end + 1 } = [ rootpath_data, 'mb0p50/cdi_rPIE_mb0p50_alpha0p000000001_oldsample_for_probe_update/' ];
N_trials_0p50( end + 1 )  = 10;

%==============
% Full Batch GD
%==============

path_data_full_gd = {};
N_trials_full_gd  = [];

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha1p00_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p10_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p01_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p001_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p00001_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p0000001_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, ];
N_trials_full_gd( end + 1 )  = 10;

path_data_full_gd{ end + 1 } = [ rootpath_data, 'full/cdi_rPIE_full_alpha0p000000001_randT_oldsample_for_probe_update/' ];
N_trials_full_gd( end + 1 )  = 10;



