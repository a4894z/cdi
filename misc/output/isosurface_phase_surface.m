function isosurface_phase_surface( obj, params )

% obj = permute( obj, [ 3, 2, 1 ] );

    szr = size( obj, 1 );
    szc = size( obj, 2 );
    szp = size( obj, 3 );

    amps    = abs( obj );
    phasesp = angle( obj ); %/phsf;%.*(amps>=0.15)/phsf; %p1a1 and p1a2

    [ x, y, z ] = meshgrid( 1 : szc, ...
                            1 : szr, ...
                            1 : szp );

    hold on;

    p = patch( isosurface( amps, params.isoslvl ));

    isonormals( x, y, z, amps, p );
    isocolors( x, y, z, phasesp, p ); %for phases
    set( p, 'FaceColor', 'interp', 'EdgeColor', 'none' )
%     set( p, 'FaceColor', 'flat', 'EdgeColor', 'none' )
    
    p2 = patch( isocaps( amps, params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
    set( p2, 'FaceColor', 'interp', 'EdgeColor', 'none' )
%     set( p2, 'FaceColor', 'flat', 'EdgeColor', 'none' )
    
%     p3 = patch( isocaps( amps( end, :, : ), params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
%     set( p3, 'FaceColor', 'interp', 'EdgeColor', 'none' )
    
    
%     xlabel('x'); xlim( [ 1, szc ])
%     ylabel('z'); ylim( [ 1, szr ])
%     zlabel('y'); zlim( [ 1, szp ])
    
    xlabel( '+x', 'fontweight', 'bold', 'fontsize', 16 ); xlim( [ 1, szc ])
    ylabel( '+z', 'fontweight', 'bold', 'fontsize', 16 ); ylim( [ 1, szr ])
    zlabel( '+y', 'fontweight', 'bold', 'fontsize', 16 ); zlim( [ 1, szp ])

    set(gca, 'xdir', 'reverse', 'fontweight','bold' );

    view( [ 45, 30 ]);

    camlight;
    lighting gouraud;
    grid on

    alpha( params.alphalvl ) 


    %     p2 = patch( isocaps( amps, params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
    %     set( p2, 'FaceColor', 'interp', 'EdgeColor', 'none' )

    hold off;
    
%     set( gca, 'xdir','reverse')

end

