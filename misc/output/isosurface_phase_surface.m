function isosurface_phase_surface( obj, params )

obj = permute( obj, [ 3, 2, 1 ] );



%     szr = size( obj, 1 );
%     szc = size( obj, 2 );
%     szp = size( obj, 3 );
%     
%     [ x, y, z ] = meshgrid( 1 : szc, ...
%                             1 : szr, ...
%                             1 : szp );
                        
                        
    amps    = abs( obj );
    phasesp = angle( obj ); %/phsf;%.*(amps>=0.15)/phsf; %p1a1 and p1a2

    if ~isfield( params, 'x_coord_length' ), params.x_coord_length = 1 : size( obj, 2 ); end
    if ~isfield( params, 'y_coord_length' ), params.y_coord_length = 1 : size( obj, 3 ); end
    if ~isfield( params, 'z_coord_length' ), params.z_coord_length = 1 : size( obj, 1 ); end
    
    if ~isfield( params, 'x_lim' ), params.x_lim = [ 1, size( obj, 2 ) ]; end
    if ~isfield( params, 'y_lim' ), params.y_lim = [ 1, size( obj, 3 ) ]; end
    if ~isfield( params, 'z_lim' ), params.z_lim = [ 1, size( obj, 1 ) ]; end
    
    [ x, y, z ] = meshgrid( params.x_coord_length, ...
                            params.z_coord_length, ...
                            params.y_coord_length );
    

                        

    hold on;

    p = patch( isosurface( x, y, z, amps, params.isoslvl ));

    isonormals( x, y, z, amps, p );
    isocolors( x, y, z, phasesp, p ); %for phases
    set( p, 'FaceColor', 'interp', 'EdgeColor', 'none' )
%     set( p, 'FaceColor', 'flat', 'EdgeColor', 'none' )
    
    p2 = patch( isocaps( x, y, z, amps, params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
    set( p2, 'FaceColor', 'interp', 'EdgeColor', 'none' )
%     set( p2, 'FaceColor', 'flat', 'EdgeColor', 'none' )
    
%     p3 = patch( isocaps( amps( end, :, : ), params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
%     set( p3, 'FaceColor', 'interp', 'EdgeColor', 'none' )
    
% 
%     xlabel( '+x', 'fontweight', 'bold', 'fontsize', 16 ); xlim( [ 1, szc ])
%     ylabel( '+z', 'fontweight', 'bold', 'fontsize', 16 ); ylim( [ 1, szr ])
%     zlabel( '+y', 'fontweight', 'bold', 'fontsize', 16 ); zlim( [ 1, szp ])
    
    xlabel( '+x', 'fontweight', 'bold', 'fontsize', 16 ); xlim( params.x_lim )
    ylabel( '+z', 'fontweight', 'bold', 'fontsize', 16 ); ylim( params.z_lim )
    zlabel( '+y', 'fontweight', 'bold', 'fontsize', 16 ); zlim( params.y_lim )
    
%     xlabel( '+x', 'fontweight', 'bold', 'fontsize', 16 ); xlim( [ min(params.x_coord_length), max(params.x_coord_length) ])
%     ylabel( '+z', 'fontweight', 'bold', 'fontsize', 16 ); ylim( [ min(params.z_coord_length), max(params.z_coord_length) ])
%     zlabel( '+y', 'fontweight', 'bold', 'fontsize', 16 ); zlim( [ min(params.y_coord_length), max(params.y_coord_length) ])

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

