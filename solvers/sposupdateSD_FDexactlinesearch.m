function [ spos_rs ] = sposupdateSD_FDexactlinesearch( sol, expt )

% figure; imagesc( abs(sol.probe.P(:,:,3)))
% 
% 
% figure; imagesc( abs(sol.phiA(:,:,3)))
% 
% figure; imagesc( abs(sol.sample.TF))
% [ tmp0 ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol );
% figure; imagesc( abs(tmp0))

%{

sol.spos.update.dpx = 5.0; % for exact line search, scan pos resolution
sol.spos.update.maxpx = 50; % max pixels to look ahead in the line search

%}


dpx         = sol.spos.update.dpx;
linesearch  = sol.spos.update.linesearch;
shifttype   = sol.spos.update.shifttype;

f_els   = zeros( length( linesearch ), length( sol.spos.rs ), 'single' );
gf      = zeros( length( sol.spos.rs ), 2, 'single' );

% x = randsample(n,10)';
% numelements = round(0.1*length(a)); 
% indices = randperm(length(a),numelements); 
% b = a(indices) 


%==================================================================================================

wNN = 0 * 45000;
wMAX = 0;
 
% sol.spos.update.indx = randperm( length( sol.spos.indxsubset ), 0.25 * length( sol.spos.indxsubset ));
% sol.spos.update.indx = sol.spos.updateorder;
 
for ss = sol.spos.update.indx
    
%     fprintf( [ num2str( [ ss, sol.spos.N ], 'Scan Position Correction %d / %d' ), '\n' ]);





    rsALL       = sol.spos.rs;
    rs0         = expt.spos.rs( ss, : );
    indxsubset  = expt.spos.indxsubset;
    
    
    
    
    
    
    spos = sol.spos.rs(  ss, : ) + [ 0, 0 ];
%     [ f0, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, 'val', shifttype );
    [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
    f0 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
    
    spos = sol.spos.rs(  ss, : ) + [ dpx, 0 ];
%     [ f1, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, 'val', shifttype );
    [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
    f1 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
    
    spos = sol.spos.rs(  ss, : ) + [ 0, dpx ];
%     [ f2, ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, 'val', shifttype );
    [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
    f2 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
    
    gfr = f1 - f0;
    gfc = f2 - f0;
    gf( ss, : ) = [ gfr, gfc  ];
    
    
    
    
    
    
    
    
    
%     spos = sol.spos.rs(  ss, : ) + 0.5 * [ dpx, 0 ];
%     [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
%     f0 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
%     
%     spos = sol.spos.rs(  ss, : ) - 0.5 * [ dpx, 0 ];
%     [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
%     f1 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
%     
%     spos = sol.spos.rs(  ss, : ) + 0.5 * [ 0 , dpx ];
%     [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
%     f2 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
%     
%     spos = sol.spos.rs(  ss, : ) - 0.5 * [ 0 , dpx ];
%     [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
%     f3 = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
%     
%     gfr = f0 - f1;
%     gfc = f2 - f3;
%     gf( ss, : ) = [ gfr, gfc  ];
    
    
    
    
    
    
    
    
%     [ ~, g ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), rs, 'grad', shifttype );
%     gf( ss, : ) = [ g( 1 ) ,g( 2 )  ];
    

    gf( ss, : ) = gf( ss, : ) / norm( gf( ss, : ), 2 );

    for ll = 1 : length( linesearch )

        spos = sol.spos.rs(  ss, : ) - linesearch( ll ) * gf( ss, : );
%         [ f_els( ll, ss ), ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, 'val', shifttype );
        [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( sol.probe.P, sol.sample.TF, sol.sample.vs, expt.meas.SI( ss ), spos, rs0, rsALL, ss, indxsubset, 'val', shifttype );
        f_els( ll, ss ) = cf.f + wMAX * cf.fMAX + wNN * cf.fNN;
    
    end
    
    6;
%{
figure
plot( linesearch, f_els( :, ss ), '-x' ); 
grid on
close all;
%}



%     [ ~, I] = min( f_els( :, ss ) );
%     sol.spos.rs( ss, : ) =  sol.spos.rs(  ss, : ) - linesearch( I ) * gf( ss, : );
    
    
% close all;
% ss

end



%==================================================================================================

% spos_rs = zeros( size( sol.spos.rs ));
spos_rs = sol.spos.rs;

if strcmp( sol.spos.update.grad, 'indiv' )
    
    % update each scan position individually
%     for ss = 1 : length( sol.spos.rs )  
    for ss = sol.spos.update.indx
        
        [ ~, I ] = min( f_els( :, ss ));
        spos_rs( ss, : ) = sol.spos.rs( ss, : ) - linesearch( I ) * gf( ss, : );

    end

elseif strcmp( sol.spos.update.grad, 'all' )
    
    % !!!! DETERMINE MINS FOR EACH INDIVIDUALLY AND TAKE AVERAGE?
    
    
    % update all scan positions using a single step length simultaneously
%     tmp0 = sum( f_els, 2 ) / length( sol.spos.rs );
    tmp0 = sum( f_els, 2 ) / length( sol.spos.update.indx );
    
    [ ~, I ] = min( tmp0 );
    aalpha = linesearch( I );

    for ss = sol.spos.update.indx
        
        spos_rs( ss, : ) = sol.spos.rs( ss, : ) - aalpha * gf( ss, : );

    end

else
    
    % 2D GRID SEARCH?
    
    
    error('!!!!!!!!!!!!!!!!!!!!!!!!!');
    
end

%==================================================================================================
% Reset scan positions if they drift too far away from guessed experimental starting points

tmp0 = ( sol.spos.update.maxcorrectr < abs( spos_rs( :, 1 ) - expt.spos.rs( :, 1 )));
tmp1 = ( sol.spos.update.maxcorrectc < abs( spos_rs( :, 2 ) - expt.spos.rs( :, 2 )));

spos_rs( tmp0, 1 ) = 0.0 * sol.spos.rs( tmp0, 1 ) + 1.0 * expt.spos.rs( tmp0, 1 );
spos_rs( tmp1, 2 ) = 0.0 * sol.spos.rs( tmp1, 2 ) + 1.0 * expt.spos.rs( tmp1, 2 );

spos_rs( tmp0, 1 ) = spos_rs( tmp0, 1 ) + 3 * ( 2 * rand( size( spos_rs( tmp0, 1 ))) - 1 );
spos_rs( tmp1, 2 ) = spos_rs( tmp1, 2 ) + 3 * ( 2 * rand( size( spos_rs( tmp1, 2 ))) - 1 );

%==================================================================================================
% Shift the scan positions so that they're on average as closest to guessed experimental starting points

% tmp0 = expt.spos.rs;
% mr = mean( spos_rs( :, 1 ) - tmp0( :, 1 ));
% mc = mean( spos_rs( :, 2 ) - tmp0( :, 2 ));
% spos_rs = spos_rs - [ mr, mc ];


