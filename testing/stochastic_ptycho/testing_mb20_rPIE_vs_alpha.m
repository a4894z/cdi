% 
%{

clear; close all; testing_mb20_rPIE_vs_alpha

%}

%====================================================================================================================================================

addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi/misc/output/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/mb0p20/';

%====================================================================================================================================================

load( [ rootpath_data, 'metrics.mat' ] );

y_lim = [-1, 6];

metrics_plot = [ 1, 12, 14, 16 ];
% metrics_plot = [ 1, 6, 14, 16 ];
% metrics_plot = [ 1, 6, 11, 12, 13, 14, 15, 16 ];

%===========================================
% using shaded error bars for trial variance
%===========================================

% skip  = 1;
% tmp0 = {};
% 
% facecolors = [ [ 0.0, 0.0, 0.0 ]; ...
%                [ 0.7, 0.0, 0.0 ]; ...
%                [ 0.0, 0.0, 0.7 ]; ...
%                [ 0.0, 0.7, 0.0 ]; ...
%                [ 0.0, 0.7, 0.7 ]; ...
%                [ 0.7, 0.0, 0.7 ]; ...
%                [ 0.7, 0.7, 0.0 ]; ...
%              ];
%          
% h1 = figure();     
% set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
%     
% 
% hold on
%     
% for jj = 1 : length( metrics_plot )
%     
%     ii = metrics_plot( jj );
%     
%     x     = metrics{ ii }.it.mtot;
%     x     = x( 1 : skip : end );
%     
% %     xconf = [ x, x( end : -1 : 1) ]; 
% %     yconf = [ metrics{ ii }.log10_meas_all_max, fliplr( metrics{ ii }.log10_meas_all_min ) ];
% % 
% % %     p = fill( xconf, yconf, 'red' );
% % %     p.FaceColor = [ 1, 0.8, 0.8 ];   
% %     
% % 
% % 
% % %     p = fill( xconf, yconf, facecolors( jj, : ) );
% %     p = fill( xconf, yconf, [ 0, 0, 0 ] );
% %     p.FaceColor = facecolors( jj, : );      
% %     p.EdgeColor = 'none';     
% %     p.FaceAlpha = 0.250;
% % %     set( h1, 'facealpha', 0.25 )
%     
% 
% 
%     plot( x, metrics{ ii }.log10_meas_all_avg( 1 : skip : end ), '-', ...
%                                                                  'linewidth', 2, ...
%                                                                  'color', [ facecolors( jj, : ), 0.75 ] )
% 
% 
%     xlabel('Epoch')
%     ylabel( { [ metrics{ ii }.name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
% 
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     
%     tmp0{ end + 1 } = metrics{ ii }.name_data; 
% 
%     %     legend( 'Location', 'northeast' ) 
% 
% end
% 
% hold off
% 
% legend( tmp0, 'interpreter', 'none' )
% 
% hold on
% for jj = 1 : length( metrics_plot )
%     
%     ii = metrics_plot( jj );
%     
%     x     = metrics{ ii }.it.mtot;
%     x     = x( 1 : skip : end );
%     xconf = [ x, x( end : -1 : 1) ]; 
%     
%     y_max = metrics{ ii }.log10_meas_all_max;
%     y_min = fliplr( metrics{ ii }.log10_meas_all_min );
%     yconf = [ y_max( 1 : skip : end ), y_min( 1 : skip : end ) ];
% 
% %     p = fill( xconf, yconf, 'red' );
% %     p.FaceColor = [ 1, 0.8, 0.8 ];   
%     
% 
% 
% %     p = fill( xconf, yconf, facecolors( jj, : ) );
%     p = fill( xconf, yconf, [ 0, 0, 0 ] );
%     p.FaceColor = facecolors( jj, : );      
%     p.EdgeColor = 'none';     
%     p.FaceAlpha = 0.08;
% %     set( h1, 'facealpha', 0.25 )
%     
% end
% 
% hold off

%==========================
% using "normal" error bars
%==========================

tmp0 = {};

skip = 5;

h1 = figure();     
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

facecolors = [ [ 0.0, 0.0, 0.0 ]; ...
               [ 0.7, 0.0, 0.0 ]; ...
               [ 0.0, 0.0, 0.7 ]; ...
               [ 0.0, 0.7, 0.0 ]; ...
               [ 0.0, 0.7, 0.7 ]; ...
               [ 0.7, 0.0, 0.7 ]; ...
               [ 0.7, 0.7, 0.0 ]; ...
             ];
         
for jj = 1 : length( metrics_plot )
    
    ii = metrics_plot( jj );
    
    x     = metrics{ ii }.it.mtot;
    x     = x( 1 : skip : end );
    
    y     = metrics{ ii }.log10_meas_all_avg;
    y     = y( 1 : skip : end );

    %========

    hold on
    
    err_max = transpose( metrics{ ii }.log10_meas_all_max( 1 : skip : end ) ) - y;

    err_min = y - transpose( metrics{ ii }.log10_meas_all_min( 1 : skip : end ) );

%     h = errorbar( x, y, err_min, err_max, '-o',            ...
%                                           'linewidth', 2, ...
%                                           'color', [ facecolors( jj, : ), 0.15 ] );
                                      
    h = errorbar( x, y, err_min, err_max, '-o',            ...
                                          'linewidth', 2, ...
                                          'color', [ facecolors( jj, : ), 0.15 ] );
                                      
    hold off
    
    %========
        
%     hold on
%     
%     S = std( metrics{ ii }.log10_meas_all, 1, 2 );
%     
%     h = errorbar( x, y, S, 'o',            ...
%                            'linewidth', 2, ...
%                            'color', [ facecolors( jj, : ), 0.15 ] );
%                      
%     hold off

    %========
    
    xlabel('Epoch')
    ylabel( { [ metrics{ ii }.name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )

    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    
    tmp0{ end + 1 } = metrics{ ii }.name_data; 

    %     legend( 'Location', 'northeast' ) 

    alpha = 0.15;
    set( [ h.Bar, h.Line ], 'ColorType', 'truecoloralpha', 'ColorData', [ h.Line.ColorData( 1 : 3 ); 255 * alpha ])

end


legend( tmp0, 'interpreter', 'none' )
tmp0 = {};

%===============================================================
% using CUSTOM FUNCTION for shaded error bars for trial variance
%===============================================================

tmp0 = {};

skip = 1;

h1 = figure();     
set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
              
facecolors = [ [ 0.0, 0.0, 0.0 ]; ...
               [ 0.8, 0.0, 0.0 ]; ...
               [ 0.0, 0.0, 0.8 ]; ...
               [ 0.0, 0.8, 0.0 ]; ...
               [ 0.5, 0.5, 0.5 ]; ...
               [ 0.8, 0.0, 0.8 ]; ...
               [ 0.8, 0.8, 0.0 ]; ...
               [ 0.0, 0.8, 0.8 ]; ...
             ];
             
for jj = 1 : length( metrics_plot )
    
    ii = metrics_plot( jj );
    
    x     = metrics{ ii }.it.mtot;
    x     = x( 1 : skip : end );
    
    y     = metrics{ ii }.log10_meas_all_avg;
    y     = y( 1 : skip : end );

    %========

    hold on
    
    err_max = transpose( metrics{ ii }.log10_meas_all_max( 1 : skip : end ) ) - y;
    err_min = y - transpose( metrics{ ii }.log10_meas_all_min( 1 : skip : end ) );

    h = shadedErrorBar( x, y, transpose( [ err_max, err_min ] ), 'lineprops', { '-', 'LineWidth', 2, 'Color', facecolors( jj, : ) }, ...
                                                                 'transparent', true, 'patchSaturation', 0.10 );     
        
                                                             
%     S = std( metrics{ ii }.log10_meas_all, 1, 2 );  
% %     S = S - y;
%     
%     err_max = y + S;
%     err_min = y - S;
%     
%     h = shadedErrorBar( x, y, transpose( [ S, S ] ), 'lineprops', { '-', 'LineWidth', 2, 'Color', facecolors( jj, : ) }, ...
%                                                                  'transparent', true, 'patchSaturation', 0.10 );                                                               
                                                             
                                                             
        
                                                             
                     
    set( h.edge, 'LineWidth', 1.75, 'LineStyle', ':' )                                                
             
%     if mod( jj, 2 )
% 
%     %     h.patch.FaceColor = 0 * [ 0.5, 0.25, 0.25 ];
% %         set( h.edge, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', facecolors( jj, : ) )
%         set( h.edge, 'LineWidth', 1.5, 'LineStyle', '--' )
%     %     h.mainLine.LineWidth = 3;
%     %     h.mainLine.LineWidth = 3;
%     else
%         
%         set( h.edge, 'LineWidth', 1.25, 'LineStyle', '-', 'Color', facecolors( jj, : ) )
%         
%     end
    
%     set( h.edge, 'transparent', true, 'patchSaturation', 0.10 )
    
    hold off
    
    %========
    
    xlabel('Epoch')
    ylabel( { [ metrics{ ii }.name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )

    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    
    tmp1 = metrics{ ii }.name_data;
    tmp1( end - 14 : end ) = [];
    tmp0{ end + 1 } = tmp1; 
    
    %     legend( 'Location', 'northeast' ) 

%     alpha = 0.15;
%     set( [ h.Bar, h.Line ], 'ColorType', 'truecoloralpha', 'ColorData', [ h.Line.ColorData( 1 : 3 ); 255 * alpha ])

end

legend( tmp0, 'interpreter', 'none' )





return

%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha1p00/independenttrials_11Aug2021_t083622/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p80/independenttrials_15Aug2021_t234638/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p60/independenttrials_15Aug2021_t061743/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p40/independenttrials_14Aug2021_t130526/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p20/independenttrials_13Aug2021_t202125/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p10/independenttrials_13Aug2021_t084624/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p09/independenttrials_12Aug2021_t224717/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p08/independenttrials_12Aug2021_t130610/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p07/independenttrials_12Aug2021_t033840/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p06/independenttrials_11Aug2021_t181332/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p05/independenttrials_10Aug2021_t071308/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p04/independenttrials_09Aug2021_t214747/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p03/independenttrials_09Aug2021_t081023/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p02/independenttrials_06Aug2021_t175923/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p01/independenttrials_04Aug2021_t145912/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p001/independenttrials_30Aug2021_t001710/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p0001/independenttrials_28Aug2021_t222641/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p00001/independenttrials_27Aug2021_t205443/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p000001/independenttrials_03Sep2021_t101824/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p0000001/independenttrials_04Sep2021_t025601/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p00000001/independenttrials_04Sep2021_t194342/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p000000001/independenttrials_05Sep2021_t103602/' ];
% N_trials( end + 1 )  = 10;

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha1p00_oldsample_for_probe_update/independenttrials_19Oct2021_t191344/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p80_oldsample_for_probe_update/independenttrials_18Oct2021_t065039/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p60_oldsample_for_probe_update/independenttrials_16Oct2021_t185322/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p40_oldsample_for_probe_update/independenttrials_15Oct2021_t065542/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p20_oldsample_for_probe_update/independenttrials_13Oct2021_t185555/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p10_oldsample_for_probe_update/independenttrials_12Oct2021_t073340/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p08_oldsample_for_probe_update/independenttrials_11Oct2021_t080320/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p06_oldsample_for_probe_update/independenttrials_10Oct2021_t112614/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p04_oldsample_for_probe_update/independenttrials_09Oct2021_t144538/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p02_oldsample_for_probe_update/independenttrials_08Oct2021_t175541/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p01_oldsample_for_probe_update/independenttrials_21Sep2021_t205815/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p001_oldsample_for_probe_update/independenttrials_21Sep2021_t205346/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p0001_oldsample_for_probe_update/independenttrials_24Sep2021_t121129/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p00001_oldsample_for_probe_update/independenttrials_24Sep2021_t224259/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p000001_oldsample_for_probe_update/independenttrials_25Sep2021_t091618/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_mb0p20_alpha0p0000001_oldsample_for_probe_update/independenttrials_25Sep2021_t195104/' ];
N_trials( end + 1 )  = 10;

%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    metrics{ jj } = load_and_plot( path_data{ jj }, N_trials( jj ) );   
    
    metrics{ jj }.rootpath_data = rootpath_data;                        

end

%====================================================================================================================================================

function metrics = load_and_plot( path_data, N_trials )

    y_lim = [-1, 6];
    
    %========
    
    sim_ptycho2DTPA = cell( N_trials, 1 );
    
    for ii = 1 : N_trials

        sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
        
    end

    %========

    meas_all_avg = 0;
    
    for ii = 1 : N_trials
        
        log10_meas_all( :, ii ) = log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); %#ok<AGROW>

        meas_all_avg = meas_all_avg + log10_meas_all( :, ii );

    end

    meas_all_avg = meas_all_avg / N_trials;

    log10_meas_all_min = transpose( min( log10_meas_all, [], 2 ));
    log10_meas_all_max = transpose( max( log10_meas_all, [], 2 ));  

    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_all( :, ii )          = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )               = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )                   = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )           = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
        metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
        
    end
    
    name_data = num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
        
    metrics.log10_meas_all_max  = log10_meas_all_max;
    metrics.log10_meas_all_min  = log10_meas_all_min;
    metrics.log10_meas_all      = log10_meas_all;
    metrics.log10_meas_all_avg  = meas_all_avg;
    metrics.name_data           = name_data;
    metrics.N_trials            = N_trials;
    metrics.path_data           = path_data;
    
    %============================================================
    % create figures of the specified metrics over the trials run
    %============================================================
    
%     skip = 1;
%     
%     h1 = figure();  
%     set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
%     x     = metrics.it( 1 ).mtot( 1 : skip : end );
%     xconf = [ x, x( end : -1 : 1) ]; 
%     yconf = [ log10_meas_all_max, fliplr( log10_meas_all_min ) ];
%     
%     p = fill( xconf, yconf, 'red' );
%     p.FaceColor = [ 1, 0.8, 0.8 ];      
%     p.EdgeColor = 'none';           
% 
%     hold on
% 
%     plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
% %     legend( 'Location', 'northeast' ) 
%     
%     hold off
     
    %========
    
    skip = 1;

    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )

    hold on

    for ii = 1 : N_trials

        plot( metrics.it( ii ).mtot( 1 : skip : end ), log10( metrics.meas_all( 1 : skip : end, ii )  ), '-', 'linewidth', 2 )

    end

    plot( metrics.it( 1 ).mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )

    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'northeast' ) 
    
    hold off

end















% %========
% 
% for jj = 1 : length( path_data )
%     
%     load_and_plot( path_data{ jj }, N_trials( jj ) )
%     
%     5;
%     
% end
% 
% %====================================================================================================================================================
% 
% function load_and_plot( path_data, N_trials )
% 
%     y_lim = [-1, 6];
%     
%     %========
%     
%     sim_ptycho2DTPA = cell( N_trials, 1 );
%     
%     for ii = 1 : N_trials
% 
%         sim_ptycho2DTPA{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
%     end
% 
%     %========
% 
%     meas_all_avg = 0;
% 
%     for ii = 1 : N_trials
% 
%         meas_all_avg = meas_all_avg + log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all ); 
% 
%     end
% 
%     meas_all_avg = meas_all_avg / N_trials;
% 
%     %========
% 
%     skip = 1;
% 
% %     figure; 
%     h1 = figure();  
%     set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
% 
%     hold on
% 
%     for ii = 1 : N_trials
% 
%         plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
%     end
% 
%     plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
%     
%     name_data = num2str( [ sim_ptycho2DTPA{ii}.sol.rPIE_alpha, sim_ptycho2DTPA{ii}.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
%     
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'northeast' ) 
% 
% 
% end
% 
% %====================================================================================================================================================
