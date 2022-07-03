%

clearvars -except expt sol

expt.spos.indxsubset = expt.spos.indx;
sol.spos.indxsubset  = expt.spos.indxsubset;

sol.spos.rs0 = sol.spos.rs;

%========================
% introduce random errors
%========================

% sol.spos.rs = gpuArray( sol.spos.rs + single( 10 * ( -1 + 2 * rand( [ sol.spos.N, 2 ] ) )));

%=======================
% SHEAR AND SCALE ERRORS
%=======================

figure; 
plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-400, 400])
ylim([-400, 400])
daspect([1 1 1])  
grid on




% % shear_x_err = sign( -2 * rand + 1 ) * 0.40;
% % shear_y_err = sign( -2 * rand + 1 ) * 0.40;
% shear_x_err = ( -2 * rand + 1 ) * 0.10;
% shear_y_err = ( -2 * rand + 1 ) * 0.10;
% 
% % scale_x_err = +1.0 + sign( -2 * rand + 1 ) * 0.40;
% % scale_y_err = +1.0 + sign( -2 * rand + 1 ) * 0.40;
% scale_x_err = +1.0 + ( -2 * rand + 1 ) * 0.10;
% scale_y_err = +1.0 + ( -2 * rand + 1 ) * 0.10;
% 
% shear_x = [ [ 1, 0 ]; [ shear_x_err, 1 ] ];
% shear_y = [ [ 1, shear_y_err ]; [ 0, 1 ] ];
% 
% scale_xy = [ [ scale_y_err, 0 ]; [ 0, scale_x_err ] ];   
% 
% sol.spos.rs = transpose( shear_x * transpose( sol.spos.rs ));            
% sol.spos.rs = transpose( shear_y * transpose( sol.spos.rs ));
% 
% sol.spos.rs = transpose( scale_xy * transpose( sol.spos.rs ));


% figure; 
% plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-400, 400])
% ylim([-400, 400])
% daspect([1 1 1])  
% grid on

%================
% ROTATION ERRORS
%================

% ???



%========

% sol.spos.rs0_err = sol.spos.rs;

%==============================================================
% zero pad sample array to account for scan position correction
%==============================================================

% expt.sample.T     = padarray( expt.sample.T, [ 0, 0 ], 1 );
% expt.sample.sz.sz = size( expt.sample.T );
% expt.sample.sz.r  = expt.sample.sz.sz( 1 );
% expt.sample.sz.c  = expt.sample.sz.sz( 2 );
% 
% expt.sample.vs.r = single( round( ( 0.5 * ( expt.sample.sz.r - expt.sz.r ) + 1 ) : ( 0.5 * ( expt.sample.sz.r + expt.sz.r ))));
% expt.sample.vs.c = single( round( ( 0.5 * ( expt.sample.sz.c - expt.sz.c ) + 1 ) : ( 0.5 * ( expt.sample.sz.c + expt.sz.c ))));

% phi = gpuArray( single( expt.probe.phi ));
% T0  = gpuArray( single( expt.sample.T  ));
% 
% spos_opt.sz      = gpuArray( single( expt.sz.sz        ));
% spos_opt.szTF    = gpuArray( single( expt.sample.sz.sz ));
% spos_opt.rc      = gpuArray( single( expt.sz.rc        ));
% spos_opt.sqrt_rc = gpuArray( single( expt.sz.sqrt_rc   ));
% spos_opt.vs_r    = gpuArray( single( expt.sample.vs.r  ));
% spos_opt.vs_c    = gpuArray( single( expt.sample.vs.c  ));

%========


sol.sample.T     = padarray( sol.sample.T, [ 0, 0 ], 1 );
sol.sample.sz.sz = size( sol.sample.T );
sol.sample.sz.r  = sol.sample.sz.sz( 1 );
sol.sample.sz.c  = sol.sample.sz.sz( 2 );

sol.sample.vs.r = single( round( ( 0.5 * ( sol.sample.sz.r - sol.sz.r ) + 1 ) : ( 0.5 * ( sol.sample.sz.r + sol.sz.r ))));
sol.sample.vs.c = single( round( ( 0.5 * ( sol.sample.sz.c - sol.sz.c ) + 1 ) : ( 0.5 * ( sol.sample.sz.c + sol.sz.c ))));


phi = gpuArray( single( sol.probe.phi ));
T0  = gpuArray( single( sol.sample.T  ));

spos_opt.sz      = gpuArray( single( sol.sz.sz        ));
spos_opt.szTF    = gpuArray( single( sol.sample.sz.sz ));
spos_opt.rc      = gpuArray( single( sol.sz.rc        ));
spos_opt.sqrt_rc = gpuArray( single( sol.sz.sqrt_rc   ));
spos_opt.vs_r    = gpuArray( single( sol.sample.vs.r  ));
spos_opt.vs_c    = gpuArray( single( sol.sample.vs.c  ));









% spos_opt.noise_model = 'poisson';
spos_opt.noise_model = 'gaussian';

    

for aa = 1 : 10
    
%     pct_of_spos_to_use   = 0.10;
%     tmp0                 = randperm( length( sol.spos.indxsubset ), round( pct_of_spos_to_use * length( sol.spos.indxsubset ))); 
%     rs_ind               = sol.spos.indxsubset( tmp0 );

    rs_ind = sol.spos.indxsubset( 1 : 1 : end );
    
    
    
    
    
    
    nDeq0   = gpuArray( not(    expt.meas.Deq0( :, :, rs_ind ) ));
    measD   = gpuArray( single( expt.meas.D( :, :, rs_ind )    ));
    if strcmp( spos_opt.noise_model, 'poisson' ), measD  = measD .^ 2; end

    spos_opt.Nspos = gpuArray( single( length( rs_ind )));

    
    rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
    
    
    [ opt_shear_y, opt_shear_x ] = scanpositions_update_2DTPA_shearxy_gridsearch( rs, T0, phi, measD, nDeq0, spos_opt );
    
    
    
    
    
    
    [ opt_scale_y, opt_scale_x ] = scanpositions_update_2DTPA_scalexy_gridsearch( rs, T0, phi, measD, nDeq0, spos_opt );

    
%     
%     rs( :, 1 ) = gpuArray( single( sol.spos.rs( rs_ind, 2 ) ));
%     rs( :, 2 ) = gpuArray( single( sol.spos.rs( rs_ind, 1 ) ));
% 
%     [ rot ] = scanpositions_update_2DTPA_rotation( rs, T0, phi, measD, nDeq0, spos_opt );
% 
% 
% %     rot = -1 * rot;
% 
%     rotT = gpuArray( [ [ cosd( rot ), -sind( rot ) ] ; ...
%                        [ sind( rot ), +cosd( rot ) ] ] );     
% 
%                
%     rs( :, 1 ) = gpuArray( single( sol.spos.rs( :, 2 ) ));
%     rs( :, 2 ) = gpuArray( single( sol.spos.rs( :, 1 ) ));           
%                
%     rs = transpose( rotT * transpose( rs ));
%     
%     sol.spos.rs( :, 2 ) = rs( :, 1 );
%     sol.spos.rs( :, 1 ) = rs( :, 2 );
    
    
    
    
    
    
    if mod( aa, 1 ) == 0

        figure; 
        plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
        set( gca, 'xdir', 'reverse' )
        set( gca, 'ydir', 'normal' )
        xlabel('xh, lab frame'); 
        ylabel('yv, lab frame');
        xlim([-400, 400])
        ylim([-400, 400])
        daspect([1 1 1])  
        grid on

        5;
    
    end

end



5;






% close all;
% 
% 
% 
% spos_opt.Naalpha_rs = gpuArray( 5 );  
% spos_opt.aalpha_rs  = gpuArray( linspace( 0.5, 5.00, spos_opt.Naalpha_rs ));   
% 
% % sol.spos.rs = sol.spos.rs + 5 * rand( size( sol.spos.rs ), 'single');
% 
% % spos_opt.optimize_rs_GD = 'individual';
% spos_opt.optimize_rs_GD = 'collective';
% 
% % spos_opt.noise_model = 'poisson';
% spos_opt.noise_model = 'gaussian';
% 
% for ii = 1 : 10
%     
%     
%     
%     
% %     pct_of_spos_to_use   = 0.1;
% %     tmp0                 = randperm( length( expt.spos.indxsubset ), round( pct_of_spos_to_use * length( expt.spos.indxsubset ))); 
% %     rs_ind               = expt.spos.indxsubset( tmp0 );
% 
%     rs_ind = expt.spos.indxsubset( 1 : 1 : end );
%     
%     nDeq0   = gpuArray( not(    expt.meas.Deq0( :, :, rs_ind ) ));
%     measD   = gpuArray( single( expt.meas.D( :, :, rs_ind )    ));
%     if strcmp( spos_opt.noise_model, 'poisson' ), measD2  = measD .^ 2; end
% 
%     spos_opt.Nspos   = gpuArray( single( length( rs_ind )));
% 
%     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
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
%     [ rs ] = scanpositions_update_2DTPA( rs, T0, phi, measD, nDeq0, spos_opt );
%     
%     
%     
%     sol.spos.rs( rs_ind, : ) = gather( rs );
%     
%     
%     
%     
%     
%     
%     if mod( ii, 1 ) == 0
%         
%         figure; 
%         plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
%         set( gca, 'xdir', 'reverse' )
%         set( gca, 'ydir', 'normal' )
%         xlabel('xh, lab frame'); 
%         ylabel('yv, lab frame');
%         xlim([-400, 400])
%         ylim([-400, 400])
%         daspect([1 1 1])  
%         grid on
%         
%         5;
% 
%     end
% 
% end






%=============================================
% exact line search for affine transformations
%=============================================

% % spos_opt.noise_model = 'poisson';
% spos_opt.noise_model = 'gaussian';
% 
% 
% % spos_opt.T = [ [ 1, 0 ]; ...
% %                [ 0, 1 ] ];
%            
% 
% spos_opt.delta_FD = 1e-2;
% 
% spos_opt.shear_x_FD = gpuArray( [ [ 1, 0 ]; [ spos_opt.delta_FD, 1 ] ] );
% spos_opt.shear_y_FD = gpuArray( [ [ 1, spos_opt.delta_FD ]; [ 0, 1 ] ] );
% spos_opt.scale_x_FD = gpuArray( [ [ 1, 0 ]; [ 0, 1 + spos_opt.delta_FD ] ] );     
% spos_opt.scale_y_FD = gpuArray( [ [ 1 + spos_opt.delta_FD, 0 ]; [ 0, 1 ] ] );             
%      
%                       
%                                 
%   
% spos_opt.affineT = [ 1 + 0.00 * ( 2 * rand - 1 ), 0.00 * ( 2 * rand - 1 ), 0.00 * ( 2 * rand - 1 ), 1 + 0.00 * ( 2 * rand - 1 ) ];
% spos_opt.affineT = gpuArray( spos_opt.affineT ); 
% 
% spos_opt.Naalpha_affineT = gpuArray( 5 );  
% spos_opt.aalpha_affineT  = gpuArray( linspace( 0.0, 0.5, spos_opt.Naalpha_affineT ));   


%=====================
% RESET SCAN POSITIONS
%=====================

% sol.spos.rs = sol.spos.rs0;
% sol.spos.rs = sol.spos.rs0_err;

% close all;
%     
% 
% figure; 
% plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
% set( gca, 'xdir', 'reverse' )
% set( gca, 'ydir', 'normal' )
% xlabel('xh, lab frame'); 
% ylabel('yv, lab frame');
% xlim([-500, 500])
% ylim([-500, 500])
% daspect([1 1 1])  
% grid on
% 
% 
% for ii = 1 : 50
%     
%     
%     
% %     pct_of_spos_to_use   = 0.10;
% %     tmp0                 = randperm( length( sol.spos.indxsubset ), round( pct_of_spos_to_use * length( sol.spos.indxsubset ))); 
% %     rs_ind               = sol.spos.indxsubset( tmp0 );
% 
%     rs_ind = sol.spos.indxsubset( 1 : 1 : end );
%     
%     nDeq0   = gpuArray( not(    expt.meas.Deq0( :, :, rs_ind ) ));
%     measD   = gpuArray( single( expt.meas.D( :, :, rs_ind )    ));
%     if strcmp( spos_opt.noise_model, 'poisson' ), measD2  = measD .^ 2; end
% 
%     spos_opt.Nspos = gpuArray( single( length( rs_ind )));
% 
%     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
% 
%     %================================================================================================================================================
% 
%     [ ~, spos_opt ] = scanpositions_update_2DTPA_affine( rs, T0, phi, measD, nDeq0, spos_opt );
% 
%     genT_aalpha = [ spos_opt.affineT( 1 ), spos_opt.affineT( 2 ); spos_opt.affineT( 3 ), spos_opt.affineT( 4 ) ]; 
%     sol.spos.rs = transpose( genT_aalpha * transpose( sol.spos.rs ));    
%     
% %     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
% 
% %     spos_opt.affineT = [ 1, 0, 0, 1 ];
%     spos_opt.affineT = [ 1 + 0.01 * sign(2 * rand - 1), 0.01 * sign(2 * rand - 1), 0.01 * sign(2 * rand - 1), 1 + 0.01 * sign(2 * rand - 1) ];
%     
%     %================================================================================================================================================
%     
%     
%     
%     
%     
%     %================================================================================================================================================
%     
%     
%     
%     
%     if mod( ii, 10 ) == 0
%         
%         figure; 
%         plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
% %         plot_2Dscan_positions( expt.spos.rs, [], rs, [] )
%         set( gca, 'xdir', 'reverse' )
%         set( gca, 'ydir', 'normal' )
%         xlabel('xh, lab frame'); 
%         ylabel('yv, lab frame');
%         xlim([-500, 500])
%         ylim([-500, 500])
%         daspect([1 1 1])  
%         grid on
%         
%         5;
% 
%     end
% 
% end




























%=====================================
% exact line search for rotation angle
%=====================================



spos_opt.delta_rot_FD = 1;
spos_opt.rot_FD       = gpuArray( [ [ sind( spos_opt.delta_rot_FD ), +cosd( spos_opt.delta_rot_FD ) ]; ...
                                    [ cosd( spos_opt.delta_rot_FD ), -sind( spos_opt.delta_rot_FD ) ] ] );     

spos_opt.Naalpha_rot = gpuArray( 15 );  
spos_opt.aalpha_rot  = gpuArray( linspace( 0.00, 0.5, spos_opt.Naalpha_rot )); 



          



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%=====================
% RESET SCAN POSITIONS
%=====================

% sol.spos.rs = sol.spos.rs0;
% sol.spos.rs = sol.spos.rs0_err;



spos_opt.shear_x = gpuArray( 0.0 );  
spos_opt.shear_y = gpuArray( 0.0 );  

spos_opt.scale_y = gpuArray( 1.0 );  
spos_opt.scale_x = gpuArray( 1.0 );  

spos_opt.Naalpha_shear = gpuArray( 10 );  
spos_opt.aalpha_shear  = gpuArray( linspace( 0.01, 0.20, spos_opt.Naalpha_shear ));        

spos_opt.Naalpha_scale = gpuArray( 10 );  
spos_opt.aalpha_scale  = gpuArray( linspace( 0.01, 0.20, spos_opt.Naalpha_scale ));   

spos_opt.Naalpha_rot = gpuArray( 10 );  
spos_opt.aalpha_rot  = gpuArray( linspace( 0.01, 0.20, spos_opt.Naalpha_rot ));   


close all;
    
figure; 
plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
set( gca, 'xdir', 'reverse' )
set( gca, 'ydir', 'normal' )
xlabel('xh, lab frame'); 
ylabel('yv, lab frame');
xlim([-400, 400])
ylim([-400, 400])
daspect([1 1 1])  
grid on


for ii = 1 : 10
    
    
    pct_of_spos_to_use   = 0.50;
    tmp0                 = randperm( length( expt.spos.indxsubset ), round( pct_of_spos_to_use * length( expt.spos.indxsubset ))); 
    rs_ind               = expt.spos.indxsubset( tmp0 );

%     rs_ind = expt.spos.indxsubset( 1 : 4 : end );
    
    nDeq0   = gpuArray( not(    expt.meas.Deq0( :, :, rs_ind ) ));
    measD   = gpuArray( single( expt.meas.D( :, :, rs_ind )    ));
    if strcmp( spos_opt.noise_model, 'poisson' ), measD2  = measD .^ 2; end

    spos_opt.Nspos = gpuArray( single( length( rs_ind )));

    rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
    
    
    
    

    

    shear_x = optimize_spos_shear_x( rs, T0, phi, measD, nDeq0, spos_opt );

    shear_x = [ [ 1,   0   ]; ...
                [ shear_x, 1 ] ];

%     rs = transpose( shear_x * transpose( rs ));
    sol.spos.rs = transpose( shear_x * transpose( sol.spos.rs ));   
    rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
    
    %========

%     shear_y = optimize_spos_shear_y( rs, T0, phi, measD, nDeq0, spos_opt );
%     
%     shear_y = [ [ 1, shear_y   ]; ...
%                 [ 0, 1 ] ];
%             
% %     rs = transpose( shear_y * transpose( rs ));
%     sol.spos.rs = transpose( shear_x * transpose( sol.spos.rs ));   
%     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
%     

    
    %========
    
%     [ scale_y ] = optimize_spos_scale_y( rs, T0, phi, measD, nDeq0, spos_opt );
%      
%     sol.spos.rs( :, 1 ) = scale_y * sol.spos.rs( :, 1 );
%     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
     
    %========
    
%     [ scale_x ] = optimize_spos_scale_x( rs, T0, phi, measD, nDeq0, spos_opt );
%      
%     sol.spos.rs( :, 2 ) = scale_x * sol.spos.rs( :, 2 );
%     rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));

    %================================================================================================================================================
    
%     rs = optimize_spos_rs( rs, T0, phi, measD, nDeq0, spos_opt );
    
    %================================================================================================================================================
    
    if mod( ii, 1 ) == 0
        
        figure; 
        plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
        set( gca, 'xdir', 'reverse' )
        set( gca, 'ydir', 'normal' )
        xlabel('xh, lab frame'); 
        ylabel('yv, lab frame');
        xlim([-400, 400])
        ylim([-400, 400])
        daspect([1 1 1])  
        grid on
        
        5;

    end

end





close all;



spos_opt.Naalpha_rs = gpuArray( 5 );  
spos_opt.aalpha_rs  = gpuArray( linspace( 0.5, 5.00, spos_opt.Naalpha_rs ));   

% sol.spos.rs = sol.spos.rs + 5 * rand( size( sol.spos.rs ), 'single');

spos_opt.optimize_rs_GD = 'individual';
% spos_opt.optimize_rs_GD = 'collective';

for ii = 1 : 10
    
    
    
    
    pct_of_spos_to_use   = 0.1;
    tmp0                 = randperm( length( expt.spos.indxsubset ), round( pct_of_spos_to_use * length( expt.spos.indxsubset ))); 
    rs_ind               = expt.spos.indxsubset( tmp0 );

%     rs_ind = expt.spos.indxsubset( 1 : 1 : end );
    
    nDeq0   = gpuArray( not(    expt.meas.Deq0( :, :, rs_ind ) ));
    measD   = gpuArray( single( expt.meas.D( :, :, rs_ind )    ));
    if strcmp( spos_opt.noise_model, 'poisson' ), measD2  = measD .^ 2; end

    spos_opt.Nspos   = gpuArray( single( length( rs_ind )));

    rs = gpuArray( single( sol.spos.rs( rs_ind, : ) ));
    


    
    
    
    
    
    
    
    
    
    [ rs ] = scanpositions_update_2DTPA( rs, T0, phi, measD, nDeq0, spos_opt );
    
    
    
    sol.spos.rs( rs_ind, : ) = gather( rs );
    
    
    
    
    
    
    if mod( ii, 1 ) == 0
        
        figure; 
        plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
        set( gca, 'xdir', 'reverse' )
        set( gca, 'ydir', 'normal' )
        xlabel('xh, lab frame'); 
        ylabel('yv, lab frame');
        xlim([-400, 400])
        ylim([-400, 400])
        daspect([1 1 1])  
        grid on
        
        5;

    end

end






    

%====================================================================================================================================================

function [ scale_y ] = optimize_spos_scale_y( rs, T0, phi, measD, nDeq0, spos_opt )

    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
    
    %========
    
    rs_scaley = transpose( spos_opt.scale_y_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_scaley = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );

    %========

    grad_Lg_scaley = sign( Lg_0 - Lg_scaley );

    %========
    
    Lg_scaley_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_scale, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_scale

        scale_y_alpha = spos_opt.scale_y + spos_opt.aalpha_scale( aa ) * grad_Lg_scaley;
        
        rs_scaley = [ scale_y_alpha * rs( :, 1 ), rs( :, 2 ) ];

        ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

        I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        sqrt_I_e = sqrt( I_e );

        Lg_scaley_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;

    end

    [ ~, II ] = min( Lg_scaley_aalpha );
    
    scale_y = spos_opt.scale_y + spos_opt.aalpha_scale( II ) * grad_Lg_scaley;
 
end


%====================================================================================================================================================

function [ scale_x ] = optimize_spos_scale_x( rs, T0, phi, measD, nDeq0, spos_opt )

    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
    
    %========
    
    rs_scalex = transpose( spos_opt.scale_x_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_scalex, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_scalex = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );

    %========

    grad_Lg_scalex = sign( Lg_0 - Lg_scalex );

    %========
    
    Lg_scalex_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_scale, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_scale

        scale_x_alpha = spos_opt.scale_x + spos_opt.aalpha_scale( aa ) * grad_Lg_scalex;
        
        rs_scalexy = [ rs( :, 1 ), scale_x_alpha * rs( :, 2 ) ];

        ind = get_indices_2Dframes( rs_scalexy, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

        I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        sqrt_I_e = sqrt( I_e );

        Lg_scalex_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;

    end

    [ ~, II ] = min( Lg_scalex_aalpha );
    
    scale_x = spos_opt.scale_x + spos_opt.aalpha_scale( II ) * grad_Lg_scalex;
 
end


%====================================================================================================================================================
% 
% function [ scale_y, scale_x ] = optimize_spos_scale_xy( rs, T0, phi, measD, nDeq0, spos_opt )
% 
%     ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     
%     %========
%     
%     rs_scalex = transpose( spos_opt.scale_x_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_scalex, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_scalex = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_scalex = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
%     
%     rs_scaley = transpose( spos_opt.scale_y_FD * transpose( rs ));
% 
%     ind = get_indices_2Dframes( rs_scaley, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%     psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%     psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%     sqrt_I_e = sqrt( I_e );
% 
%     Lg_scaley = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% %     Lp_scaley = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
% 
%     %========
% 
%     grad_Lg_scalex = Lg_0 - Lg_scalex;
%     grad_Lg_scaley = Lg_0 - Lg_scaley;
%     grad_Lg        = [ grad_Lg_scaley, grad_Lg_scalex ];
%     grad_Lg_N      = grad_Lg ./ sqrt( sum( abs( grad_Lg ).^ 2, 2 ));
%     
%     %========
%     
%     Lg_scalexy_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha, 'single' ));
%     
%     for aa = 1 : spos_opt.Naalpha
% 
%         scale_xy_alpha = [ spos_opt.scale_y, spos_opt.scale_x ] + spos_opt.aalpha( aa ) * grad_Lg_N;
%         
%         scale_xy  = [ [ scale_xy_alpha( 1 ), 0                   ]; ...
%                       [ 0,                   scale_xy_alpha( 2 ) ] ];
%   
%         rs_scalexy = transpose( scale_xy * transpose( rs ));
% 
%         ind = get_indices_2Dframes( rs_scalexy, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
%         psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
%         psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
% 
%         I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
%         sqrt_I_e = sqrt( I_e );
% 
%         Lg_scalexy_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;
% 
%     end
% 
%     [ ~, II ] = min( Lg_scalexy_aalpha );
%     
%     scale_xy = [ spos_opt.scale_y, spos_opt.scale_x ] + spos_opt.aalpha( II ) * grad_Lg_N;
%     
%     scale_y = scale_xy( 1 );
%     scale_x = scale_xy( 2 );
%     
% end

%====================================================================================================================================================

function shear_x = optimize_spos_shear_x( rs, T0, phi, measD, nDeq0, spos_opt )

    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
    
    %========
    
    rs_shearx = transpose( spos_opt.shear_x_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_shearx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_shearx = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_shearx = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );

    %========

    grad_Lg_shearx = sign( Lg_0 - Lg_shearx );

    %========
    
    Lg_shearx_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_shear, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_shear

        s_x_alpha = spos_opt.shear_x + spos_opt.aalpha_shear( aa ) * grad_Lg_shearx;
        
        shear_x  = [ [ 1,          0 ]; ...
                     [ s_x_alpha,  1 ] ];
  
        rs_shearx = transpose( shear_x * transpose( rs ));

        ind = get_indices_2Dframes( rs_shearx, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

        I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        sqrt_I_e = sqrt( I_e );

        Lg_shearx_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;

    end

    [ ~, II ] = min( Lg_shearx_aalpha );
    
    shear_x = spos_opt.shear_x + spos_opt.aalpha_shear(  II ) * grad_Lg_shearx;

end


%====================================================================================================================================================

function shear_y = optimize_spos_shear_y( rs, T0, phi, measD, nDeq0, spos_opt )

    ind = get_indices_2Dframes( rs, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_0 = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_0 = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
    
    %========
    
    rs_sheary = transpose( spos_opt.shear_y_FD * transpose( rs ));

    ind = get_indices_2Dframes( rs_sheary, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
    psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

    I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
    sqrt_I_e = sqrt( I_e );

    Lg_sheary = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );
%     Lp_shearx = sum( sum( sum( nDeq0 .* ( I_e - measD2 .* log( I_e )), 1 ), 2 )) / ( spos_opt.Nspos * spos_opt.rc );

    %========

    grad_Lg_sheary = sign( Lg_0 - Lg_sheary );

    %========
    
    Lg_sheary_aalpha = gpuArray( zeros( 1, spos_opt.Naalpha_shear, 'single' ));
    
    for aa = 1 : spos_opt.Naalpha_shear

        s_y_alpha = spos_opt.shear_y + spos_opt.aalpha_shear( aa ) * grad_Lg_sheary;
        
        shear_y  = [ [ 1,  s_y_alpha ]; ...
                     [ 0,  1 ] ];
  
        rs_sheary = transpose( shear_y * transpose( rs ));

        ind = get_indices_2Dframes( rs_sheary, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
        psi = reshape( T0( ind ), [ spos_opt.sz, 1, spos_opt.Nspos ]) .* phi;
        psi = fft( fft( fftshift( fftshift( psi, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;

        I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );
        sqrt_I_e = sqrt( I_e );

        Lg_sheary_aalpha( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e ) .^ 2, 1 ), 2 ) / spos_opt.rc ) / spos_opt.Nspos;

    end

    [ ~, II ] = min( Lg_sheary_aalpha );
    
    shear_y = spos_opt.shear_y + spos_opt.aalpha_shear(  II ) * grad_Lg_sheary;

end





