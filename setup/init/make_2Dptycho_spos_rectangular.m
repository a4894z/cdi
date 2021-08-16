function [ spos_rs ] = make_2Dptycho_spos_rectangular( spos )

    % make_2DTPA_ptychospos_rectangular
    % make_ptychospos_2DTPA_rectangular
    % make_ptychospos_3DtomoTPA_rectangular
    % make_ptychospos_3DlaminoTPA_rectangular
    % make_ptychospos_3DlaminoTMS_rectangular

    %==================================================================================================

    ppR = spos.shftpx_r * (( 1 : 1 : spos.Nrow ) - 1 );    % vector of row positions 
    ppC = spos.shftpx_c * (( 1 : 1 : spos.Ncol ) - 1 );    % vector of col positions

    spos_rs = transpose( [ repelem( ppR, spos.Ncol ); repmat( ppC, 1, spos.Nrow )]);

    % give the assumed sample frame scan positions an overall translation so that the 1st position is at the origin:
    spos_rs = spos_rs - min( spos_rs, [], 1 );
    spos_rs = spos_rs - max( spos_rs, [], 1 ) * 0.5;  
    
    %==================================================================================================
    
    if ( spos.Nrow * spos.Ncol ) > 1

        tmp0 = 1 : length( spos_rs );
        tmp0 = transpose( reshape( tmp0, [ spos.Ncol, spos.Nrow ] ));
    %     tmp0 = flipud( tmp0 );

        for ii = 1 : size( tmp0, 1 )

            if ~mod( ii, 2 )

                tmp0( ii, : ) = fliplr( tmp0( ii, : ));

            end
        end

        spos_sio = transpose( tmp0 );      
        spos_sio = spos_sio( : );

        spos_rs( :, 2 ) = spos_rs( spos_sio, 2 );

    end

    %=======================
    
    % give the assumed sample frame scan positions an overall translation so that the 1st position is at the origin:
    spos_rs = spos_rs - min( spos_rs, [], 1 );
    spos_rs = spos_rs - max( spos_rs, [], 1 ) * 0.5;  
    
    %=======================
    
%     spos.indxsubset = 1 : spos.N;
%     spos_rs      = spos_rs;

    %=======================
    
    spos_rs  = single( spos_rs );
    
    %==================================================================================================
    %-------------- introduce misc probe position goofiness and scan position errors ------------------
    %==================================================================================================

    startrow = -1 * 0 + 0.0 * rand; 
    startcol = -1 * 0 + 0.0 * rand; 

    % the "ezperimenta" scan positions, without errors
    spos_rs( :, 1 ) = spos_rs( :, 1 ) + startrow;
    spos_rs( :, 2 ) = spos_rs( :, 2 ) + startcol;

    %=======================
    % translation
   
    % spos.translate( 1 ) = 6;
    % spos.translate( 2 ) = -20;
    % spos_rs( :, 1 ) = spos_rs( :, 1 ) + spos.translate( 1 );
    % spos_rs( :, 2 ) = spos_rs( :, 2 ) + spos.translate( 2 );

    %=======================
    % % shear in vertical direction
    
    % spos.shear.r = +0.3;
    % spos_rs = spos_rs * [ 1, 0; spos.shear.r, 1 ];

    % % shear in horizontal direction
    % spos.shear.c = +0.05;
    % spos_rs = spos_rs * [ 1, spos.shear.c; 0, 1 ];

    %=======================
    % % rotate
    % 
    % spos.rot = 30 * pi / 180;
    % spos_rs = spos_rs * [ cos( spos.rot ), sin( spos.rot ); -sin( spos.rot ), cos( spos.rot ) ];

    %=======================
    % scale

    % sx = 0.9;
    % sy = 1.0;
    % spos_rs = spos_rs * [ sx, 0; 0, sy ];

    %=======================
    % randomize positions

    % spos_rs = spos_rs + 0 * round( 2 * ( 2 * rand( spos.N, 2 ) - 1 ));

%     spos_rs( :, 1 ) = spos_rs( :, 1 ) + round( spos.randR * ( 2 * rand( spos.N, 1 ) - 1 ));
%     spos_rs( :, 2 ) = spos_rs( :, 2 ) + round( spos.randC * ( 2 * rand( spos.N, 1 ) - 1 ));

    %=======================

    % figure; 
    % plot_2Dscan_positions( spos_rs, ( 1 : length( spos_rs )), spos_rs, [] );
    % title(num2str( 1e6 * [ expt.csys.z2.dLy, expt.csys.z2.dLx ],'vertical pixel size = %.3e um, horizontal pixel size = %.3e um' )) 
    % %axis tight
    % 
    % 5;
    



end

%==================================================================================================
