%

%{

clear;
load( '~/Documents/MATLAB/Code/psycho_rat/run/experimental/Star_fine_mid_out.mat' );
spos_test = 517;

% testing_scanpos_finitediff_exactls;

%}

sol.errmetric.compute = 'val';

dpx = 0.5;
lnsrch = transpose( 0 : 0.5 : 3 );
f_els = zeros( length( lnsrch ), length( sol.spos.rs ), 'single' );
gf = zeros( length( sol.spos.rs ), 2, 'single' );

for ss = 1 : length( sol.spos.rs )  
    
%     fprintf( [ num2str( [ ss, sol.spos.N ], 'Scan Position Correction %d / %d' ), '\n' ]);
        
    rs = sol.spos.rs( ss, : );

    spos = rs + [ 0, 0 ];
    [ f0, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( ss ), spos, sol );

    spos = rs + [ dpx, 0 ];
    [ f1, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( ss ), spos, sol );

    spos = rs + [ 0, dpx ];
    [ f2, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( ss ), spos, sol );


    gfr = f1 - f0;
    gfc = f2 - f0;

    gf( ss, : ) = [ gfr, gfc  ];
    gf( ss, : ) = gf( ss, : ) / norm( gf( ss, : ), 2 );

    for ll = 1 : length( lnsrch )

        spos = rs - lnsrch( ll ) * gf( ss, : );
        [ f_els( ll, ss ), ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( ss ), spos, sol );

    end
    
%{
figure
plot( lnsrch, f_els( :, ss ), '-x' ); 
grid on
close all;
%}



%     [ ~, I] = min( f_els( :, ss ) );
%     sol.spos.rs( ss, : ) = rs - lnsrch( I ) * gf( ss, : );
    
    
% close all;
% ss

end

%==================================================================================================

if strcmp( sol.spos.update.grad, 'indiv' )
    
    % update each scan position individually
    for ss = 1 : length( sol.spos.rs )  

        [ ~, I ] = min( f_els( :, ss ));
        sol.spos.rs( ss, : ) = sol.spos.rs( ss, : ) - lnsrch( I ) * gf( ss, : );

    end

elseif strcmp( sol.spos.update.grad, 'all' )
    
    % !!!! DETERMINE MINS FOR EACH INDIVIDUALLY AND TAKE AVERAGE?
    
    
    % update all scan positions using a single step length simultaneously
    tmp0 = sum( f_els, 2 ) / length( sol.spos.rs );

    [ ~, I ] = min( tmp0 );
    aalpha = lnsrch( I );

    for ss = 1 : length( sol.spos.rs )  

        sol.spos.rs( ss, : ) = sol.spos.rs( ss, : ) - aalpha * gf( ss, : );

    end

end

%==================================================================================================
% SHIFT THE SPOS SO THAT THEY'RE ON AVERAGE AS CLOSEST TO STARTING SPOS AS POSSIBLE

tmp0 = expt.spos.rs( expt.spos.indxsubset, : );
mr = mean( sol.spos.rs( :, 1 ) - tmp0( :, 1 ));
mc = mean( sol.spos.rs( :, 2 ) - tmp0( :, 2 ));
sol.spos.rs = sol.spos.rs - [ mr, mc ];

%==================================================================================================
% RESET SPOS IF THEY STRAY TOO FAR AWAY FROM STARTING POINT

tmp0 = expt.spos.rs( expt.spos.indxsubset, : );
Ireset = ( 10 < sqrt(( sol.spos.rs( :, 1 ) - tmp0( :, 1 )).^2 + ( sol.spos.rs( :, 2 ) - tmp0( :, 2 )).^2 ));
% sol.spos.rs( Ireset == 1, : ) = tmp0( Ireset == 1, : );
sol.spos.rs( Ireset, : ) = tmp0( Ireset, : );






%{

clear; close all; load /home/ash/Desktop/18oct_9x9_r1/Star_fine_mid_out_1000.mat;


figure
plot_2Dscan_positions( double( expt.spos.rs( expt.spos.indxsubset, : )), ...
                       [], ...
                       double( sol.spos.rs ), ...
                       expt.spos.indxsubset );
                 

tmp0 = expt.spos.rs( expt.spos.indxsubset, : );               
mr = mean( sol.spos.rs( :, 1 ) - tmp0( :, 1 ));
mc = mean( sol.spos.rs( :, 2 ) - tmp0( :, 2 ));
                               
sol.spos.rs = sol.spos.rs - [ mr, mc ];

mr = mean( sol.spos.rs( :, 1 ) - tmp0( :, 1 ))
mc = mean( sol.spos.rs( :, 2 ) - tmp0( :, 2 ))
%}



% figure
% plot( lnsrch, f_els, '-' ); 
% hold on
% plot( lnsrch, tmp0, '-x','linewidth', 2 ); 
% hold off
% 
% grid on
% 
% 5;






