function ptycho2DTPA_mkimg_sample_SCPM( sol, expt )


        close all;
        
        h1 = figure();  
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        plot_2Dscan_positions( expt.spos.rs, [], sol.spos.rs, [] )
%         plot_2Dscan_positions( expt.spos.rs, [], rs, [] )
        set( gca, 'xdir', 'reverse' )
        set( gca, 'ydir', 'normal' )
        xlabel('xh, lab frame'); 
        ylabel('yv, lab frame');
%         xlim([-500, 500])
%         ylim([-500, 500])
        daspect([1 1 1])  
        grid on
        
        export_fig( num2str( sol.it.epoch, 'scanpositions_%d.jpg' ), '-r120.0' )
        close all;
    
        %==============================
        % probe mode correlation matrix
        %==============================

%         P = reshape( sol.probe.phi, [ sol.sz.rc, sol.probe.scpm.N ] );
% 
%         corr_matrix_scpm = ctranspose( P ) * P;
%         
%         h1 = figure();  
%         set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
%     
%         imagesc( log10( 1 + abs( corr_matrix_scpm )))
%         axis square
%         colorbar
%         colormap( expt.cm.blj )
%         title('Probe Mode Correlation Matrix')
%         
%         export_fig( num2str( sol.it.epoch, 'corr_matrix_scpm_%d.jpg' ), '-r120.0' )
%         close all;
        
        %================
        % sample plotting
        %================

%         hsv_img = complex2hsv( sol.sample.T, 0, 0, 0 );
%         rgb_img = hsv2rg( hsv_img );
        
        close all;
    
        pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sample.sz.c);
        pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
        pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
        
        pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sample.sz.r);
        pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
        pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
        
        h1 = figure();  
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        ax1 = subplot(131);
        imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.sample.T ), [ 0, 1.00 ]); 

%         daspect([1 1 1]); 
        axis square

        colorbar
        colormap( ax1, expt.cm.blj )
        grid on; 
        title('abs sample')
        

        ax2 = subplot(132);
        imagesc( pltopts.xaxis, pltopts.yaxis, angle( sol.sample.T ), [ -pi, pi ] ); 
        
%         daspect([1 1 1]); 
        axis square

        colorbar
        colormap( ax2, expt.cm.hsvD )
        grid on; 
        title('phase sample')
        
        
        subplot(133)
        imagescHSV( sol.sample.T, pltopts  ); 

%         daspect([1 1 1]); 
        axis square

        grid on;
        title('HSV ( V = mag, H = phs ) sample')
        
        export_fig( num2str( sol.it.epoch, 'sample_%d.jpg' ), '-r120.0' )
        close all;
     
        %==============
        % SCPM plotting
        %==============

        [ scpm ] = compute_scpm_photonocc( sol.probe.phi );

        close all;
        pltopts.xaxis = expt.csys.z2.dLx * (1 : sol.sz.c);
        pltopts.xaxis = pltopts.xaxis - min( pltopts.xaxis );
        pltopts.xaxis = pltopts.xaxis - 0.5 * max( pltopts.xaxis );
        
        pltopts.yaxis = expt.csys.z2.dLy * (1 : sol.sz.r);
        pltopts.yaxis = pltopts.yaxis - min( pltopts.yaxis );
        pltopts.yaxis = pltopts.yaxis - 0.5 * max( pltopts.yaxis );
        
        h1 = figure();        
        set( h1, 'Visible', 'off', 'Position', [ 1, 1, 1920, 1080 ] )
        
        for pp = 1 : sol.probe.scpm.N

%             subaxis(2,sol.probe.scpm.N,pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
            subplot( 2, double( sol.probe.scpm.N ), double( pp ) )   
            imagescHSV( sol.probe.phi( :, :, pp ), pltopts ); 
%             imagescHSV( log10( 1 + 15^-1 * abs( sol.probe.phi( :, :, pp ))) .* exp( 1i * angle( sol.probe.phi( :, :, pp ))), pltopts); 
            %axis square
            daspect([1 1 1]); 
            %axis off
%             title( { num2str( scpm.occ( pp ), 'HSV ( V = mag, H = phs ), occupancy = %.4f' ), num2str(  scpm.fro2TOT, 'fro2TOT = %.4e' ) })
            title( { num2str( scpm.occ( pp ), 'occupancy = %.4f' ), num2str(  scpm.fro2TOT, 'fro2TOT = %.4e' ) })
            grid on;
            set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
      
        end
        
        for pp = 1 : sol.probe.scpm.N

%             subaxis(2,sol.probe.scpm.N,sol.probe.scpm.N +pp,'SpacingVert',0,'MR',0.01, 'ML',0.01,'MT',0.18, 'MB',0.18); 
            ax( pp ) = subplot( 2, double( sol.probe.scpm.N ), double( sol.probe.scpm.N + pp ) );
            imagesc( pltopts.xaxis, pltopts.yaxis, abs( sol.probe.phi( :, :, pp ))); 
%             imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 15^-1 * abs( sol.probe.phi( :, :, pp ))) ); 
            %axis square
            daspect([1 1 1]);
%             colormap gray; 
            colormap( ax( pp ), expt.cm.blj )
            colorbar
            grid on;
            set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
            title('abs probe')
        end

%         export_fig( num2str( sol.it.exwv, 'probe_%d.jpg' ), '-r90.0' )
        export_fig( num2str( sol.it.epoch, 'probe_%d.jpg' ), '-r90.0' )
        close all;

        %=======================================================
        % probe intensity plotting, example measurement plotting
        %=======================================================   
        
        absPmodes2 = sqrt( sum( abs( sol.probe.phi ) .^ 2, 3 ));
        
        [ phi, ~ ] = enforce_2DTPAsposview( sol.probe.phi, sol.sample.T, sol.sample.vs.r, sol.sample.vs.c, sol.spos.rs( round( 0.5 * sol.spos.N ), : ), sol.spos.shifttype );
%         phi = fftshift( fft2( fftshift( phi ))) / sqrt( numel( phi ));
        V = fft( fftshift( fft( fftshift( phi, 1 ), [], 1 ), 2 ), [], 2 ) / sqrt( numel( phi ));
        V = fftshift( sqrt( sum( abs( V ) .^ 2, 3 )));

        h1 = figure();        
        set( h1, 'Visible', 'off', 'Position',[ 1, 1, 1920, 1080 ] )
        
        a1 = subplot(121);
        imagesc( pltopts.xaxis, pltopts.yaxis, absPmodes2 ); 
%         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
        daspect([1 1 1]);
        colormap( a1, expt.cm.blj ); 
        colorbar
        grid on;
        set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('abs probe')

        a2 = subplot(122);
        imagesc( log10( 1 + 10^1 * abs( V )));
%         imagesc( pltopts.xaxis, pltopts.yaxis, log10( 1 + 1e-0 * absPmodes2 )); 
        daspect([1 1 1]);
        colormap( a2, expt.cm.blj ); 
        colorbar
        grid on;
        set( gca, 'GridColor', [0.8, 0.0, 0.0], 'GridLineStyle', '--', 'GridAlpha', 0.5 )
        title('abs^2 fft2 of typical exit wave')
        
        
%         export_fig( num2str( sol.it.exwv, 'absPmodes2_%d.jpg' ), '-r90.0' )
        export_fig( num2str( sol.it.epoch, 'absPmodes2_%d.jpg' ), '-r90.0' )
        close all;
        

end

