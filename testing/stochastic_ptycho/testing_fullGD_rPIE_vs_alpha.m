% 

%
%{

clear; close all; testing_fullGD_rPIE_vs_alpha

%}

%====================================================================================================================================================

% addpath( '/net/s8iddata/export/8-id-ECA/Analysis/atripath/shadedErrorBar/' );  

%====================================================================================================================================================

rootpath_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/minibatch_vs_stoch_vs_full_rPIE_alpha/full/';

%====================================================================================================================================================

%%{

load( [ rootpath_data, 'metrics.mat' ] );

y_lim = [-1, 6];

% metrics_plot = [ 1, 2, 3, 4, 5, 6, 7 ];
% metrics_plot = [ 1, 2, 5, 6 ];
metrics_plot = [ 1, 4, 5, 6 ];

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
% using "normal "error bars
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
%     tmp1( end - 14 : end ) = [];
    tmp0{ end + 1 } = tmp1; 
    
    %     legend( 'Location', 'northeast' ) 

%     alpha = 0.15;
%     set( [ h.Bar, h.Line ], 'ColorType', 'truecoloralpha', 'ColorData', [ h.Line.ColorData( 1 : 3 ); 255 * alpha ])

end

legend( tmp0, 'interpreter', 'none' )



return

%}

%====================================================================================================================================================

path_data = {};
N_trials  = [];

%========

% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p000000001_randT/independenttrials_04Sep2021_t141134/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p00000001_randT/independenttrials_02Sep2021_t105630/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p0000001_randT_a/independenttrials_31Aug2021_t092025/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p0000001_randT_b/independenttrials_01Sep2021_t192318/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p000001_randT_b/independenttrials_01Sep2021_t035527/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p000001_randT_a/independenttrials_31Aug2021_t092152/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p00001_randT/independenttrials_31Aug2021_t091637/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p0001_randT/independenttrials_26Aug2021_t033321/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p001_randT/independenttrials_25Aug2021_t115703/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p005_randT/independenttrials_26Aug2021_t190514/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p01_randT/independenttrials_14Aug2021_t081939/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p02_randT/independenttrials_17Aug2021_t085244/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p03_randT/independenttrials_18Aug2021_t015318/' ];
% N_trials( end + 1 )  = 18;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p04_randT/independenttrials_19Aug2021_t085857/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p05_randT/independenttrials_15Aug2021_t004151/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p06_randT/independenttrials_20Aug2021_t021023/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p07_randT/independenttrials_20Aug2021_t191435/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p08_randT/independenttrials_21Aug2021_t123327/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p09_randT/independenttrials_22Aug2021_t054834/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p10_randT/independenttrials_15Aug2021_t170501/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p20_randT/independenttrials_18Aug2021_t095541/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p30_randT/independenttrials_16Aug2021_t092814/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p40_randT/independenttrials_19Aug2021_t031108/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p50_randT/independenttrials_17Aug2021_t015208/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p60_randT/independenttrials_19Aug2021_t202608/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p70_randT/independenttrials_20Aug2021_t133027/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p80_randT/independenttrials_21Aug2021_t064520/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p90_randT/independenttrials_22Aug2021_t000627/' ];
% N_trials( end + 1 )  = 10;
% 
% path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha1p00_randT/independenttrials_13Aug2021_t155848/' ];
% N_trials( end + 1 )  = 10;

%========

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha1p00_randT_oldsample_for_probe_update/independenttrials_21Oct2021_t022747/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p10_randT_oldsample_for_probe_update/independenttrials_02Oct2021_t235705/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p01_randT_oldsample_for_probe_update/independenttrials_20Sep2021_t142149/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p001_randT_oldsample_for_probe_update/independenttrials_02Oct2021_t071559/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p00001_randT_oldsample_for_probe_update/independenttrials_01Oct2021_t143516/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p0000001_randT_oldsample_for_probe_update/independenttrials_30Sep2021_t215726/' ];
N_trials( end + 1 )  = 10;

path_data{ end + 1 } = [ rootpath_data, '/cdi_rPIE_full_alpha0p000000001_randT_oldsample_for_probe_update/independenttrials_30Sep2021_t051649/' ];
N_trials( end + 1 )  = 10;


%====================================================================================================================================================

for jj = 1 : length( path_data )
    
    metrics{ jj } = load_and_plot( path_data{ jj }, N_trials( jj ) );   %#ok<SAGROW>
    
    metrics{ jj }.rootpath_data = rootpath_data;                        %#ok<SAGROW>
    
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

%     name_data = num2str( [ sim_ptycho2DTPA{ 1 }.sol.rPIE_alpha, sim_ptycho2DTPA{ 1 }.sol.spos.rand_spos_subset_pct ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');

    %=================================================
    % create a struct from metrics data for future use
    %=================================================
    
    for ii = 1 : N_trials
        
        metrics.meas_all( :, ii )          = sim_ptycho2DTPA{ ii }.sol.metrics.meas_all;
        metrics.timing( ii )               = sim_ptycho2DTPA{ ii }.sol.timings;
        metrics.it( ii )                   = sim_ptycho2DTPA{ ii }.sol.it;
        metrics.rPIE_alpha( ii )           = sim_ptycho2DTPA{ ii }.sol.rPIE_alpha;
%         metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
        
    end
    
%     name_data = num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
    name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');
    
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








%{

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
%         metrics.rand_spos_subset_pct( ii ) = sim_ptycho2DTPA{ ii }.sol.spos.rand_spos_subset_pct;
        
    end
    
    name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');
    
    metrics.meas_all_avg_of_log = meas_all_avg;
    metrics.name_data           = name_data;
    
    %=====================
    % create metrics plots
    %=====================

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
% %     name_data = num2str( [ metrics.rPIE_alpha( 1 ), metrics.rand_spos_subset_pct( 1 ) ], 'rPIE_alpha = %0.8f, MBpct = %0.4f');
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

%     figure; 
    h1 = figure();  
    set( h1, 'Visible', 'on', 'Position',[ 1, 1, 1920, 1080 ] )
        
    hold on

    for ii = 1 : N_trials

        plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), log10( sim_ptycho2DTPA{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )

    end

    plot( sim_ptycho2DTPA{ ii }.sol.it.mtot( 1 : skip : end ), meas_all_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
       
    xlabel('Epoch')
    ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
    hold off
    title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
    grid on
    ylim( y_lim )
    legend( 'Location', 'southwest' ) 


end

%}

%====================================================================================================================================================



















% 
% function load_and_plot( path_data, N_trials )
% 
%     y_lim = [ -1, 6 ];
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
%     name_data = num2str( sim_ptycho2DTPA{ii}.sol.rPIE_alpha, 'rPIE_alpha = %0.8f');
%         
%     xlabel('Epoch')
%     ylabel( { [ name_data, ',' ], 'Cost Function Value' }, 'Interpreter', 'none' )
%     hold off
%     title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
%     grid on
%     ylim( y_lim )
%     legend( 'Location', 'southwest' ) 
% 
% 
% end
% 
% %====================================================================================================================================================
% 
% 
% 



























% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% y_lim = [-1, 6];
% 
% %======================================================================================================================
% % using random sample starts, full batch update order doesn't matter so we need to introduce stochasticity in some form
% %======================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha100_randT/independenttrials_13Aug2021_t155848/';
%     cdi_ER_rPIE_0p100{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p100_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p100_avg = cdi_ER_rPIE_0p100_avg + log10( cdi_ER_rPIE_0p100{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p100_avg = cdi_ER_rPIE_0p100_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p100{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p100{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
%     
% end
% 
% plot( cdi_ER_rPIE_0p100{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p100_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 1.00, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% legend( 'Location', 'southwest' )
% 
% ylim( y_lim )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p100')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha90_randT/independenttrials_22Aug2021_t000627/';
%     cdi_ER_rPIE_0p90{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p90_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p90_avg = cdi_ER_rPIE_0p90_avg + log10( cdi_ER_rPIE_0p90{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p90_avg = cdi_ER_rPIE_0p90_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p90{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p90{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p90{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p90_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.90, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p90')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha80_randT/independenttrials_21Aug2021_t064520/';
%     cdi_ER_rPIE_0p80{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p80_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p80_avg = cdi_ER_rPIE_0p80_avg + log10( cdi_ER_rPIE_0p80{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p80_avg = cdi_ER_rPIE_0p80_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p80{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p80{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p80{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p80_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.80, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p80')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha70_randT/independenttrials_20Aug2021_t133027/';
%     cdi_ER_rPIE_0p70{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p70_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p70_avg = cdi_ER_rPIE_0p70_avg + log10( cdi_ER_rPIE_0p70{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p70_avg = cdi_ER_rPIE_0p70_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p70{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p70{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p70{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p70_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.70, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p70')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha60_randT/independenttrials_19Aug2021_t202608/';
%     cdi_ER_rPIE_0p60{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p60_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p60_avg = cdi_ER_rPIE_0p60_avg + log10( cdi_ER_rPIE_0p60{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p60_avg = cdi_ER_rPIE_0p60_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p60{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p60{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p60{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p60_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.60, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p60')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha50_randT/independenttrials_17Aug2021_t015208/';
%     cdi_ER_rPIE_0p50{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p50_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p50_avg = cdi_ER_rPIE_0p50_avg + log10( cdi_ER_rPIE_0p50{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p50_avg = cdi_ER_rPIE_0p50_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p50{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p50{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p50{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p50_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.50, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p50')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha40_randT/independenttrials_19Aug2021_t031108/';
%     cdi_ER_rPIE_0p40{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p40_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p40_avg = cdi_ER_rPIE_0p40_avg + log10( cdi_ER_rPIE_0p40{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p40_avg = cdi_ER_rPIE_0p40_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p40{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p40{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p40{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p40_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.40, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p40')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha30_randT/independenttrials_16Aug2021_t092814/';
%     cdi_ER_rPIE_0p30{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p30_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p30_avg = cdi_ER_rPIE_0p30_avg + log10( cdi_ER_rPIE_0p30{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p30_avg = cdi_ER_rPIE_0p30_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p30{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p30{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p30{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p30_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.30, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p30')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha20_randT/independenttrials_18Aug2021_t095541/';
%     cdi_ER_rPIE_0p20{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p20_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p20_avg = cdi_ER_rPIE_0p20_avg + log10( cdi_ER_rPIE_0p20{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p20_avg = cdi_ER_rPIE_0p20_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p20{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p20{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p20{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p20_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.20, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p20')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha10_randT/independenttrials_15Aug2021_t170501/';
%     cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p10_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p10_avg = cdi_ER_rPIE_0p10_avg + log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p10_avg = cdi_ER_rPIE_0p10_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p10_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.10, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p10')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha09_randT/independenttrials_22Aug2021_t054834/';
%     cdi_ER_rPIE_0p09{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p09_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p09_avg = cdi_ER_rPIE_0p09_avg + log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p09_avg = cdi_ER_rPIE_0p09_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p09_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.09, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p09')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha08_randT/independenttrials_21Aug2021_t123327/';
%     cdi_ER_rPIE_0p08{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p08_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p08_avg = cdi_ER_rPIE_0p08_avg + log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p08_avg = cdi_ER_rPIE_0p08_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p08_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.08, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p08')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha07_randT/independenttrials_20Aug2021_t191435/';
%     cdi_ER_rPIE_0p07{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p07_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p07_avg = cdi_ER_rPIE_0p07_avg + log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p07_avg = cdi_ER_rPIE_0p07_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p07_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.07, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p07')
% 
% %====================================================================================================================================================
% 
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha06_randT/independenttrials_20Aug2021_t021023/';
%     cdi_ER_rPIE_0p06{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p06_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p06_avg = cdi_ER_rPIE_0p06_avg + log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p06_avg = cdi_ER_rPIE_0p06_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p06_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.06, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p06')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha05_randT/independenttrials_15Aug2021_t004151/';
%     cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p05_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p05_avg = cdi_ER_rPIE_0p05_avg + log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p05_avg = cdi_ER_rPIE_0p05_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p05_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.05, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p05')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha04_randT/independenttrials_19Aug2021_t085857/';
%     cdi_ER_rPIE_0p04{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p04_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p04_avg = cdi_ER_rPIE_0p04_avg + log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p04_avg = cdi_ER_rPIE_0p04_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p04_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.04, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p04')
% 
% %====================================================================================================================================================
% 
% % N_trials = 18;
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha03_randT/independenttrials_18Aug2021_t015318/';
%     cdi_ER_rPIE_0p03{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p03_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p03_avg = cdi_ER_rPIE_0p03_avg + log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p03_avg = cdi_ER_rPIE_0p03_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p03_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.03, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p03')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha02_randT/independenttrials_17Aug2021_t085244/';
%     cdi_ER_rPIE_0p02{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p02_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p02_avg = cdi_ER_rPIE_0p02_avg + log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p02_avg = cdi_ER_rPIE_0p02_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p02_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.02, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p02')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha01_randT/independenttrials_14Aug2021_t081939/';
%     cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p01_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg + log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p01_avg = cdi_ER_rPIE_0p01_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p01_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.01, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p01')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha0p001_randT/independenttrials_25Aug2021_t115703/';
%     cdi_ER_rPIE_0p001{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p001_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p001_avg = cdi_ER_rPIE_0p001_avg + log10( cdi_ER_rPIE_0p001{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p001_avg = cdi_ER_rPIE_0p001_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p001{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p001{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p001{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p001_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.001, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p001')
% 
% %====================================================================================================================================================
% 
% N_trials = 10;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_full_alpha0p0001_randT/independenttrials_26Aug2021_t033321/';
%     cdi_ER_rPIE_0p0001{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% %========
% 
% cdi_ER_rPIE_0p0001_avg = 0;
% 
% for ii = 1 : N_trials
%     
%     cdi_ER_rPIE_0p0001_avg = cdi_ER_rPIE_0p0001_avg + log10( cdi_ER_rPIE_0p0001{ ii }.sol.metrics.meas_all );
% 
% end
% 
% cdi_ER_rPIE_0p0001_avg = cdi_ER_rPIE_0p0001_avg / N_trials;
% 
% %========
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p0001{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p0001{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% plot( cdi_ER_rPIE_0p0001{ ii }.sol.it.mtot( 1 : skip : end ), cdi_ER_rPIE_0p0001_avg( 1 : skip : end ), '--', 'linewidth', 4, 'color', [ 0.0, 0.0, 0.0 ] )
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.0001, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %========
% 
% clear('cdi_ER_rPIE_0p0001')
% 
% 
% 
% 
% 
% 
% 
% 
% return
% 
% %=======================================================
% % just one trial for same random starts for sample/probe
% %=======================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha01_gpu4/independenttrials_04Aug2021_t150149/';
%     cdi_ER_rPIE_0p01{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p01{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p01{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.01, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha02_gpu4/independenttrials_07Aug2021_t011459/';
%     cdi_ER_rPIE_0p02{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p02{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p02{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.02, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha03_gpu4/independenttrials_07Aug2021_t174121/';
%     cdi_ER_rPIE_0p03{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p03{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p03{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.03, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha04_gpu4/independenttrials_08Aug2021_t100521/';
%     cdi_ER_rPIE_0p04{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p04{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p04{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.04, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha05_gpu4/independenttrials_09Aug2021_t023601/';
%     cdi_ER_rPIE_0p05{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p05{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p05{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.05, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha06/independenttrials_09Aug2021_t133912/';
%     cdi_ER_rPIE_0p06{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p06{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p06{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.06, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha07/independenttrials_09Aug2021_t151739/';
%     cdi_ER_rPIE_0p07{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p07{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p07{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.07, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha08/independenttrials_09Aug2021_t165553/';
%     cdi_ER_rPIE_0p08{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p08{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p08{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.08, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha09/independenttrials_09Aug2021_t183320/';
%     cdi_ER_rPIE_0p09{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p09{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p09{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.09, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% N_trials = 1;
% 
% for ii = 1 : N_trials
%     
%     path_data = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/04Aug_minibatch_vs_stoch_vs_full_rPIE_alpha/cdi_rPIE_totalGD_alpha10/independenttrials_09Aug2021_t201045/';
%     cdi_ER_rPIE_0p10{ ii } = load( [ path_data, num2str( ii, 'trial_%d/sim_ptycho2DTPA.mat') ] );
% 
% end
% 
% skip = 1;
% 
% figure; 
% 
% hold on
% 
% for ii = 1 : N_trials
%     
%     plot( cdi_ER_rPIE_0p10{ ii }.sol.it.mtot( 1 : skip : end ), log10( cdi_ER_rPIE_0p10{ ii }.sol.metrics.meas_all( 1 : skip : end ) ), '-', 'linewidth', 2 )
% 
% end
% 
% xlabel('Epoch')
% ylabel('rPIE alpha = 0.10, Cost Function Value')
% hold off
% title('$log_{10}\bigg[ \frac{1}{N_s} \sum_{s=1}^{N_s} \left \Vert \sqrt{W_s} -  \sqrt{ \sum_p \left\vert \mathcal{F}[ \phi_p \odot T_s ] \right\vert^2} \right\Vert^2_F \bigg]$', 'FontWeight','bold', 'FontSize', 14, 'Interpreter', 'latex' );
% grid on
% ylim( y_lim )
% legend( 'Location', 'southwest' )
% 
% %====================================================================================================================================================
% 
% 

