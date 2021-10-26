%% 

%
%{

clear; close all; testing_update_order_mb_vs_alpha

%}

%====================================================================================================================================================

%===========
% full batch
%===========

%==================
% blind_block_stoch
%==================

%===============
% mini-batch 50%
%===============

[ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_50pct;

%===============
% mini-batch 33%
%===============

% [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_33pct;

%===============
% mini-batch 20%
%===============

% [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_20pct;

%===============
% mini-batch 10%
%===============

%==============
% mini-batch 5%
%==============


%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    load_and_plot( path_data{ jj }, N_trials( jj ), ylabel_data{ jj } )

end

%====================================================================================================================================================

function load_and_plot( path_data, N_trials, ylabel_data )

    y_lim = [-1, 6];
    
    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );

    end

    %========

    meas_all_avg = 0;

    for ii = 1 : N_trials

        meas_all_avg = meas_all_avg + log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); 

    end

    meas_all_avg = meas_all_avg / N_trials;

    %========

    skip = 1;

%     figure; 
    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

    hold on

    for ii = 1 : N_trials

        plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

    end

    plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
    
%     name_data = num2str( [ sim_ptycho2DTPA{ii}.sol.rPIE_alpha, sim_ptycho2DTPA{ii}.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
    name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');

    name_data = [ ylabel_data, ', ', name_data ];
    
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
%     legend( 'Location', 'northeast' ) 
    legend( 'Location', 'southwest' ) 
    
end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_50pct

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p50/';

%========
% 0.00001

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p00001/independenttrials_05Sep2021_t225627/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p00001_oldsample_for_probe_update/independenttrials_26Sep2021_t065757/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p00001_recompute_exitwaves_after_sample/independenttrials_29Sep2021_t134918/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.0001

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p0001/independenttrials_05Sep2021_t054812/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p0001_oldsample_for_probe_update/independenttrials_25Sep2021_t164343/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p0001_recompute_exitwaves_after_sample/independenttrials_28Sep2021_t125138/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.001

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p001/independenttrials_03Sep2021_t171323/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p001_oldsample_for_probe_update/independenttrials_25Sep2021_t023411/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p001_recompute_exitwaves_after_sample/independenttrials_27Sep2021_t170128/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.01

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p01/independenttrials_04Sep2021_t163059/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p01_oldsample_for_probe_update/independenttrials_24Sep2021_t122324/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p01_recompute_exitwaves_after_sample/independenttrials_26Sep2021_t211212/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.02

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p02/independenttrials_08Sep2021_t114039/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p02_oldsample_for_probe_update/independenttrials_15Oct2021_t223951/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p02_recompute_exitwaves_after_sample/independenttrials_15Oct2021_t071648/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.04

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p04/independenttrials_09Sep2021_t045047/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p04_oldsample_for_probe_update/independenttrials_14Oct2021_t214205/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p04_recompute_exitwaves_after_sample/independenttrials_14Oct2021_t082019/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.06

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p06/independenttrials_09Sep2021_t215940/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p06_oldsample_for_probe_update/independenttrials_13Oct2021_t214225/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p06_recompute_exitwaves_after_sample/independenttrials_13Oct2021_t083513/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.08

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p08/independenttrials_10Sep2021_t150917/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p08_oldsample_for_probe_update/independenttrials_12Oct2021_t225140/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p08_recompute_exitwaves_after_sample/independenttrials_12Oct2021_t080224/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.10

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p10/independenttrials_07Sep2021_t163956/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p10_oldsample_for_probe_update/independenttrials_12Oct2021_t005707/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p10_recompute_exitwaves_after_sample/independenttrials_11Oct2021_t150307/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.20

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p20/independenttrials_11Sep2021_t081609/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p20_oldsample_for_probe_update/independenttrials_21Oct2021_t104045/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p50_alpha0p20_recompute_exitwaves_after_sample/independenttrials_22Oct2021_t065104/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';


end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_33pct

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p33/';

%========
% 1.00

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha1p00/independenttrials_13Aug2021_t153855/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha1p00_oldsample_for_probe_update/independenttrials_07Oct2021_t204113/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha1p00_recompute_exitwaves_after_sample/independenttrials_08Oct2021_t054423/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.80

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p80/independenttrials_18Aug2021_t225524/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p80_oldsample_for_probe_update/independenttrials_07Oct2021_t002735/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p80_recompute_exitwaves_after_sample/independenttrials_07Oct2021_t093558/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.60

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p60/independenttrials_19Sep2021_t134338/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p60_oldsample_for_probe_update/independenttrials_06Oct2021_t040251/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p60_recompute_exitwaves_after_sample/independenttrials_06Oct2021_t130432/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.40

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p40/independenttrials_19Sep2021_t043435/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p40_oldsample_for_probe_update/independenttrials_05Oct2021_t075837/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p40_recompute_exitwaves_after_sample/independenttrials_05Oct2021_t170034/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.20

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p20/independenttrials_18Sep2021_t192438/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p20_oldsample_for_probe_update/independenttrials_04Oct2021_t115500/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p20_recompute_exitwaves_after_sample/independenttrials_04Oct2021_t205727/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.10

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p10/independenttrials_16Aug2021_t140514/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p10_oldsample_for_probe_update/independenttrials_03Oct2021_t155218/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p10_recompute_exitwaves_after_sample/independenttrials_04Oct2021_t005401/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.01

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p01/independenttrials_14Aug2021_t090209/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p01_oldsample_for_probe_update/independenttrials_25Sep2021_t132343/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p01_recompute_exitwaves_after_sample/independenttrials_26Sep2021_t222157/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.001

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p001_a/independenttrials_15Aug2021_t024044/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p001_b/independenttrials_30Aug2021_t092104/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p001_oldsample_for_probe_update/independenttrials_26Sep2021_t002144/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p001_recompute_exitwaves_after_sample/independenttrials_27Sep2021_t115149/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.0001

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p0001/independenttrials_29Aug2021_t071708/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p0001_oldsample_for_probe_update/independenttrials_26Sep2021_t112108/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p33_alpha0p0001_recompute_exitwaves_after_sample/independenttrials_28Sep2021_t012323/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_20pct

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p20/';

%========
% 0.001

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p001/independenttrials_30Aug2021_t001710/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p001_oldsample_for_probe_update/independenttrials_21Sep2021_t205346/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p001_recompute_exitwaves_after_sample/independenttrials_22Sep2021_t084158/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.01

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p01/independenttrials_04Aug2021_t145912/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p01_oldsample_for_probe_update/independenttrials_21Sep2021_t205815/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p01_recompute_exitwaves_after_sample/independenttrials_22Sep2021_t084203/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========
% 0.05

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05/independenttrials_10Aug2021_t071308/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05_oldsample_for_probe_update/independenttrials_15Sep2021_t120331/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05_recompute_exitwaves_after_sample/independenttrials_15Sep2021_t120733/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';


end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_10pct

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p10/';

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05/independenttrials_10Aug2021_t195103/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05_oldsample_for_probe_update/independenttrials_15Sep2021_t114753/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p10_alpha0p05_recompute_exitwaves_after_sample/independenttrials_15Sep2021_t193044/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = minibatch_05pct

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p05/';

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01/independenttrials_16Aug2021_t173355/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01_oldsample_for_probe_update/independenttrials_24Sep2021_t123331/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p05_alpha0p01_recompute_exitwaves_after_sample/independenttrials_24Sep2021_t230058/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = fullbatch_block_stoch

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/blind_block_stoch/';

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p05/independenttrials_06Aug2021_t093611/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p05_oldsample_for_probe_update/independenttrials_16Sep2021_t151721/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p05_recompute_exitwaves_after_sample/independenttrials_17Sep2021_t075247/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p01/independenttrials_04Aug2021_t150026/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p01_oldsample_for_probe_update/independenttrials_22Sep2021_t104744/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: only Tk';
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_stochGD_alpha0p01_recompute_exitwaves_after_sample/independenttrials_23Sep2021_t012217/' ];
% N_trials( end + 1 )  = 10;
% ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

end

%====================================================================================================================================================

function [ rootpath_data, path_data, N_trials, ylabel_data ] = fullbatch_total_grad

path_data = {};
N_trials  = [];
ylabel_data = {};

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/full/';

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p05_randT/independenttrials_15Aug2021_t004151/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: mixed Tk and Tk+1';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p05_randT_oldsample_for_probe_update/independenttrials_16Sep2021_t162457/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: only Tk';

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p05_randT_recompute_exitwaves_after_sample/independenttrials_17Sep2021_t090613/' ];
N_trials( end + 1 )  = 10;
ylabel_data{ end + 1 } = 'probe update: recompute exitwaves using Tk+1';

end