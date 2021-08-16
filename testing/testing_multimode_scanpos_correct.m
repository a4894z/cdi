%

%{

clear; close all; 
restoredefaultpath
expt.path.code = '~/Documents/MATLAB/Code/cdi';   
cd( expt.path.code );
addpath(genpath(pwd));
clearvars -except expt

clear; close all; setup_forwardmodel_2DTPA;



clear; close all; testing_multimode_scanpos_correct;


delete(gcp('nocreate'))
parpool(10);


% brute force error metric landscape, and gradients using built-in matlab function
%
% sum_q | | F[ phi_s ] | - sqrt( D_s ) |^2
% vs
% sum_s sum_q | | F[ phi_s ] | - sqrt( D_s ) |^2



%}

%==================================================================================================

% new random seed
rng shuffle;

% load data
% load( '~/Documents/MATLAB/Code/cdi/run/misctesting2DTPA/simPRexpt_2DTPA.mat' );

load( '~/Documents/MATLAB/Code/cdi/Star_fine_mid_out.mat' );

%==================================================================================================

tmp0 = make_gaussian( sol.sz.sz, [ 0.5 * sol.sz.r + 1, 0.5 * sol.sz.c + 1 ], [ 0.9 * sol.sz.r, 0.9 * sol.sz.c ]);
sol.measLPF = fftshift( tmp0 );

sol.spos.updateorder = randperm( length( sol.spos.N ));
% sol.spos.updateorder = 1 : length( sol.spos.N );

% figure; 
% plot_2Dscan_positions( expt.spos.rs, ( 1 : length( expt.spos.rs )), [],  [] );

figure
% plot_2Dscan_positions( expt.spos.pyx_rbv_all(1 : 2 : length(expt.spos.indxsubset)), 1 : 2 : length(expt.spos.indxsubset), [], [] );
plot_2Dscan_positions( expt.spos.pyx_rbv_all, 1 : 5 : length( expt.spos.indxsubset ), [], [] );
grid on; 



%==================================================================================================

% w = 0.7;
% sol.probe.P = w * sol.probe.P + ( 1 - w ) * expt.probe.P;
% w = 0.7;
% sol.sample.TF = w * sol.sample.TF + ( 1 - w ) * expt.sample.TF;

% sol.spos.shifttype = 'px';
sol.spos.shifttype = 'subpx';  

% [ sol.phi ] = ERupdate_exwv( sol, expt );
% [ sol.sample.TF ] = ePIEupdate_sample( sol.probe.P, sol.sample.TF, sol.phi, sol );
% [ sol.probe.P ] = ePIEupdate_probemodes( sol.probe.P, sol.sample.TF, sol.phi, sol );
% % [ sol.probe.P ] = DMupdate_probemodes( sol.sample.TF, sol.phi, sol );

% for ss = sol.spos.updateorder
%     
%     [ phi, TFview ] = enforce_2DTPAsposview( sol.probe.P, sol.sample.TF, sol.sample.vs, sol.spos.rs(ss,:), sol.spos.shifttype );
% 
% end

%==================================================================================================

sol2.smpl_trans = sol.sample.TF;
sol2.probe = sol.probe.P;
alg.modmetrc_gam = 0.5;
N.vsy = ( 0.5 * sol.sample.sz.r - 0.5 * sol.sz.r + 1 ) : ( 0.5 * sol.sample.sz.r + 0.5 * sol.sz.r );
N.vsx = ( 0.5 * sol.sample.sz.c - 0.5 * sol.sz.c + 1 ) : ( 0.5 * sol.sample.sz.c + 0.5 * sol.sz.c );
N.sqrt_NxNy = sol.sz.sqrt_rc;
N.Nyo = sol.sample.sz.r;
N.Nxo = sol.sample.sz.c;

%==================================================================================================

% sposerr = round( 5 * ( 2 * [ rand( sol.spos.N, 1), rand( sol.spos.N, 1) ] - 1 ));
% sposerr = [ -3.2, 1.6 ];
% sol.spos.rs = sposerr + expt.spos.rs;
% sol.spos.rs = round( sposerr + expt.spos.rs );

sol.spos.shifttype = 'px';
% sol.spos.shifttype = 'subpx';  

dpxr = 1/1;
dpxc = 1/1;
spos_search_r = -30 : dpxr : 30;
spos_search_c = -30 : dpxc : 30;

fA      = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
grA     = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gcA     = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );

gr_bi   = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gc_bi   = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );

f_fd    = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gr_fd   = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gc_fd   = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );

fold    = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gfrold  = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );
gfcold  = zeros( length( spos_search_r ), length( spos_search_c ), 'single' );










% spos_study = round( sol.spos.N / 2 );
spos_study = 2571;

P = sol.probe.P;
TF = sol.sample.TF;
vs = sol.sample.vs;


% f       = gpuArray( zeros( length( spos_search_r ), length( spos_search_c ), 'single' ));
% grA     = gpuArray( zeros( length( spos_search_r ), length( spos_search_c ), 'single' ));
% gcA     = gpuArray( zeros( length( spos_search_r ), length( spos_search_c ), 'single' ));
% gc_bi   = gpuArray( zeros( length( spos_search_r ), length( spos_search_c ), 'single' ));
% gr_bi   = gpuArray( zeros( length( spos_search_r ), length( spos_search_c ), 'single' ));
% tmp1 = gpuArray( sol.probe.P );
% tmp2 = gpuArray( sol.sample.TF );

kkA = 1;
kkB = 1;



% for ss = 1 : sol.spos.N
for ss = spos_study 

    rs = sol.spos.rs( ss, : );
    meas = expt.meas.SI( ss );
%     tmp3 = gpuArray( expt.meas.SI( ss ));
    
    for ii = 1 : length( spos_search_r )
        
        tmp0 = spos_search_r( ii );
 
        for jj = 1 : length( spos_search_c )

            spos = rs + [ tmp0 , spos_search_c( jj ) ];

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            tic
            
            [ f0, gf0 ] = errormetric_2DTPAmeas_scanposgrad( P, TF, vs, meas, spos, 'val_grad', sol.spos.shifttype );

 
            fA( ii, jj ) = f0;
            grA( ii, jj ) = gf0( 1 );
            gcA( ii, jj ) = gf0( 2 );
            
            t = toc;
            tA( kkA ) = t;
            kkA = kkA + 1;
            
            
            
            tic
            
            [ f0, ~ ] = errormetric_2DTPAmeas_scanposgrad( P, TF, vs, meas, spos + [ 0, 0 ], 'val', sol.spos.shifttype );
            [ f1, ~ ] = errormetric_2DTPAmeas_scanposgrad( P, TF, vs, meas, spos + [ dpxr, 0 ], 'val', sol.spos.shifttype );
            [ f2, ~ ] = errormetric_2DTPAmeas_scanposgrad( P, TF, vs, meas, spos + [ 0, dpxc ], 'val', sol.spos.shifttype );
            
            f_fd( ii, jj ) = f0;
            gr_fd( ii, jj ) = ( f1 - f0 ) / dpxr;
            gc_fd( ii, jj ) = ( f2 - f0 ) / dpxc;
            
            t = toc;

            t_fd( kkB ) = t;
            kkB = kkB + 1;
   
            
            
            
            
            
            
            
            
            
            
            
            
            
%             rsMAX = 30; 
%             sposNN = [];
            
%             [r,c] = find( expt.spos.indxscanorder == spos_study );
%             
% %             sposNN( end+1 ) = expt.spos.indxscanorder( r, c );
%             if c + 1 < 81, sposNN( end+1 ) = expt.spos.indxscanorder( r, c + 1 ); end
%             if c - 1 > 1,  sposNN( end+1 ) = expt.spos.indxscanorder( r, c - 1 ); end
%             if r + 1 < 81, sposNN( end+1 ) = expt.spos.indxscanorder( r + 1, c ); end
%             if r - 1 > 1,  sposNN( end+1 ) = expt.spos.indxscanorder( r - 1, c ); end
%             if ( r - 1 > 1 ) && (c - 1 > 1),   sposNN( end+1 ) = expt.spos.indxscanorder( r - 1, c - 1 ); end
%             if ( r - 1 > 1 ) && (c + 1 < 81),  sposNN( end+1 ) = expt.spos.indxscanorder( r - 1, c + 1 ); end
%             if ( r + 1 < 81 ) && (c - 1 > 1),  sposNN( end+1 ) = expt.spos.indxscanorder( r + 1, c - 1 ); end
%             if ( r + 1 < 81 ) && (c + 1 < 81), sposNN( end+1 ) = expt.spos.indxscanorder( r + 1, c + 1 ); end
            
%             tmp9 = 0;
%             for aa = 1 : length( sposNN )
%                 
%                 rsNN( aa, : ) = sol.spos.rs( sposNN( aa ), : );
%                 
%                 tmp9 = tmp9 + 1 / ( 1e-4 + sum( abs( spos - rsNN( aa, : ) ).^2 ));
%                 
%             end




            indxsubset = expt.spos.indxsubset;
            cindx = ss;
            rs0 = rs;
            rsALL = sol.spos.rs;

%             tmp9 = 0;
%             for aa = 1 : length( indxsubset )
%                 
%                 if indxsubset( aa ) ~= cindx
%                     
%                     rsNN = rsALL( indxsubset( aa ), : );
%                     tmp9 = tmp9 + 1 / ( 1e-4 + sum( abs( spos - rsNN ).^1 ));
%                     
%                 end
%                 
%             end
%             
%             fNN( ii, jj ) = tmp9;
%             fMAX( ii, jj ) = sum( abs( spos - rs0 ).^2 );
            
            

            
            [ cf, ~ ] = errormetric_2DTPAmeas_scanposgrad_TESTING( P, TF, vs, meas, spos, rs0, rsALL, cindx, indxsubset, 'val', sol.spos.shifttype );
            
            fNN( ii, jj ) = cf.fNN;
            fMAX( ii, jj ) = cf.fMAX;
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
%             [ f0, gf0 ] = objf_modulus_scanposgrad( spos, sol2, expt.meas.SI( ss ).D, N, alg, expt.meas.SI( ss ).Deq0, 'both' );
%             fold( ii, jj ) = f0;
%             gfrold( ii, jj ) = gf0( 1 );
%             gfcold( ii, jj ) = gf0( 2 );





        end
        
        
        
        
        

        [ ss, ii ]

    end
    
    [ gc_bi( :, : ), gr_bi( :, : ) ] = gradient( fA( :, : ), 1 );
    
end

figure; 
plot(tA); 
hold on; 
plot(t_fd); 
hold off

figure; 
plot( tA ./ t_fd );

%==================================================================================================

% sz = sqrt( abs( gcA ) .^ 2 + abs( grA ).^2 ); 
% grA_uv = 1 * grA ./ ( 1e-7 + sz );
% gcA_uv = 1 * gcA ./ ( 1e-7 + sz );
% 
% skp = 1;
% 
% 
% 
% 
% for ss = spos_study  
% % for ss = 1 : sol.spos.N
%     
% %     sindx = sol.spos.indxnew( ss );
% %     rs = sol.spos.rs( sindx, : );
% %     
% 
%     figure;
%     imagesc( spos_search_c, spos_search_r, fA( :, : ))
%     colormap jet
%     grid on
%     title('using analytic computation of cost function')
% 
%     hold on
%     
%     plot( spos_search_c(c),spos_search_r(r), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color',[1,0,0] )
%     
%     quiver( spos_search_c(1 : skp : end), ...
%             spos_search_r(1 : skp : end), ...
%             +gcA_uv( 1 : skp : end, 1 : skp : end ), ...
%             +grA_uv( 1 : skp : end, 1 : skp : end ), ...
%             0.3, 'color', [0 0 0]); daspect([1 1 1]);
%         
% 
%     
%     hold off
% 
%     set(gca, 'ydir','normal')
%         
% end
% 
% % figure; 
% % imagesc( spos_search_c, spos_search_r, f( :, :, 1 ) + f( :, :, 2 ) + f( :, :, 3 ) + f( :, :, 4 ) )
% % colormap jet

%==================================================================================================

sz = sqrt( abs( gc_fd ) .^ 2 + abs( gr_fd ).^2 ); 
gr_fd_uv = 1 * gr_fd ./ ( 1e-7 + sz );
gc_fd_uv = 1 * gc_fd ./ ( 1e-7 + sz );

skp = 1;

[r,c]=find(fA == min(fA(:)));

for ss = spos_study  
% for ss = 1 : sol.spos.N
    
%     sindx = sol.spos.indxnew( ss );
%     rs = sol.spos.rs( sindx, : );
%     

    figure;
    imagesc( spos_search_c, spos_search_r, fA( :, : ))
    colormap jet
    grid on
    title('using FINITE DIFFERENCE computation of cost function')

    hold on
    
    plot( spos_search_c(c),spos_search_r(r), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color',[1,0,0] )
    
    quiver( spos_search_c(1 : skp : end), ...
            spos_search_r(1 : skp : end), ...
            +gc_fd_uv( 1 : skp : end, 1 : skp : end ), ...
            +gr_fd_uv( 1 : skp : end, 1 : skp : end ), ...
            0.3, 'color', [0 0 0]); daspect([1 1 1]);

    hold off

    set(gca, 'ydir','normal')
    
end

% figure; 
% imagesc( spos_search_c, spos_search_r, f( :, :, 1 ) + f( :, :, 2 ) + f( :, :, 3 ) + f( :, :, 4 ) )
% colormap jet


%==================================================================================================

wMAX = 20;
wNN = 0;

tmp0 = fA + wMAX * fMAX + wNN * fNN;
[ gc_biNN, gr_biNN ] = gradient( tmp0, 1 );

sz = sqrt( abs( gc_biNN ) .^ 2 + abs( gr_biNN ).^2 ); 
gr_biNNU = 1 * gr_biNN ./ ( 1e-7 + sz );
gc_biNNU = 1 * gc_biNN ./ ( 1e-7 + sz );

skp = 1;

[r,c]=find(tmp0 == min(tmp0(:)));

for ss = spos_study  
% for ss = 1 : sol.spos.N
    
%     sindx = sol.spos.indxnew( ss );
%     rs = sol.spos.rs( sindx, : );
%     

    figure;
    imagesc( spos_search_c, spos_search_r, tmp0 )
    colormap jet
    grid on
    title('using FINITE DIFFERENCE computation of cost function')

    hold on
    
    plot( spos_search_c(c),spos_search_r(r), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color',[1,0,0] )
    
    quiver( spos_search_c(1 : skp : end), ...
            spos_search_r(1 : skp : end), ...
            +gc_biNNU( 1 : skp : end, 1 : skp : end ), ...
            +gr_biNNU( 1 : skp : end, 1 : skp : end ), ...
            0.3, 'color', [0 0 0]); daspect([1 1 1]);

    hold off

    set(gca, 'ydir','normal')
    
end

%==================================================================================================

wMAX = 20;
wNN = 0;

tmp0 = fA + wMAX * fMAX + wNN * fNN;
[ gc_biNN, gr_biNN ] = gradient( tmp0, 1 );

sz = sqrt( abs( gc_biNN ) .^ 2 + abs( gr_biNN ).^2 ); 
gr_biNNU = 1 * gr_biNN ./ ( 1e-7 + sz );
gc_biNNU = 1 * gc_biNN ./ ( 1e-7 + sz );

skp = 1;

[r,c]=find(tmp0 == min(tmp0(:)));

for ss = spos_study  
% for ss = 1 : sol.spos.N
    
%     sindx = sol.spos.indxnew( ss );
%     rs = sol.spos.rs( sindx, : );
%     

    figure;
    imagesc( spos_search_c, spos_search_r, tmp0 )
    colormap jet
    grid on
    title('using FINITE DIFFERENCE computation of cost function')

    hold on
    
    plot( spos_search_c(c),spos_search_r(r), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color',[1,0,0] )
    
    quiver( spos_search_c(1 : skp : end), ...
            spos_search_r(1 : skp : end), ...
            +gc_biNNU( 1 : skp : end, 1 : skp : end ), ...
            +gr_biNNU( 1 : skp : end, 1 : skp : end ), ...
            0.3, 'color', [0 0 0]); daspect([1 1 1]);

    hold off

    set(gca, 'ydir','normal')
    
end


return

%==================================================================================================

% sz = sqrt( abs( gfcold ) .^ 2 + abs( gfrold ).^2 ); 
% g_rUold = 1 * gfrold ./ ( 1e-7 + sz );
% g_cUold = 1 * gfcold ./ ( 1e-7 + sz );
% 
% skp = 1;
% 
% for ss = spos_study  
% % for ss = 1 : sol.spos.N
%     
% %     sindx = sol.spos.indxnew( ss );
% %     rs = sol.spos.rs( sindx, : );
% %     
% 
%     figure;
%     imagesc( spos_search_c, spos_search_r, fold( :, : ))
%     colormap jet
%     title('using analytic computation of OLD cost function')
% 
%     hold on
% 
%     quiver( spos_search_c(1 : skp : end), ...
%             spos_search_r(1 : skp : end), ...
%             +g_cUold( 1 : skp : end, 1 : skp : end ), ...
%             +g_rUold( 1 : skp : end, 1 : skp : end ), ...
%             0.3, 'color', [0 0 0]); daspect([1 1 1]);
% 
%     hold off
%     set(gca, 'ydir','normal')
%     
% end

% figure; 
% imagesc( spos_search_c, spos_search_r, f( :, :, 1 ) + f( :, :, 2 ) + f( :, :, 3 ) + f( :, :, 4 ) )
% colormap jet

%==================================================================================================

v

% gr_bi = diff( f, 1, 1 );
% gc_bi = diff( f, 1, 2 );
% [ gc_bi, gr_bi ] = gradient( f( :, :, rsplt ), 1 );

sz = sqrt( abs( gc_bi ) .^ 2 + abs( gr_bi ).^2 ); 
gr_biU = 1 * gr_bi ./ ( 1e-7 + sz );
gc_biU = 1 * gc_bi ./ ( 1e-7 + sz );

skp = 1;

[r,c]=find(fA == min(fA(:)));

figure;
imagesc( spos_search_c, spos_search_r, fA( :, : ) )
colormap jet
grid on
title('using BUILT IN gradient() function')

hold on

plot( spos_search_c(c),spos_search_r(r), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color',[1,0,0] )
    
quiver( spos_search_c(1 : skp : end), ...
        spos_search_r(1 : skp : end), ...
        +gc_biU( 1 : skp : end, 1 : skp : end ), ...
        +gr_biU( 1 : skp : end, 1 : skp : end ), ...
        0.3, 'color', [0 0 0]); daspect([1 1 1]);

hold off
set(gca, 'ydir','normal')
    

return























%==================================================================================================

% LINESEARCH USING SINGLE POSITIONS INDEPENDENTLY VS ALL-IN-ONE-GO

% lnsrch = linspace( 0, 3, 41 );
lnsrch = 0 : 0.1 : 3;

fls = zeros( length( lnsrch ), sol.spos.N, 'single' );
    
for ss = find( sol.spos.indxnew == spos_test )
    
    sindx = sol.spos.indxnew( ss );
    rs = sol.spos.rs( sindx, : );
    
%     spos = sol.spos.rs( ss, : );

    [ f0, gf0 ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( sindx ), rs, sol );

    sz = sqrt( abs( gf0( 1 ) ) .^ 2 + abs( gf0( 2 ) ).^2 ); 
    gf0( 1 ) = 1 * gf0( 1 ) ./ ( 1e-7 + sz );
    gf0( 2 ) = 1 * gf0( 2 ) ./ ( 1e-7 + sz );



    for ll = 1 : length( lnsrch )

        spos = rs - lnsrch( ll ) * gf0;

        [ fls( ll, ss ), ~ ] = errormetric_2DTPAmeas_scanposgrad( sol.probe.P, sol.sample.TF, expt.meas.SI( sindx ), spos, sol );

    end

end


flsavg = 0;

for ii = 1 : size( fls, 2 )
    
    flsavg = flsavg + fls( :, ii );
    
end


figure; 
plot( lnsrch, fls, '-x' ); 
hold on
plot( lnsrch, flsavg / size( fls, 2 ), '-o', 'linewidth', 2 ); 
hold off
grid on
legend

5;

%==================================================================================================




































%==================================================================================================

% brute force error metric landscape, and gradients using finite differences

%==================================================================================================

% brute force error metric landscape, and gradients using analytic computation





