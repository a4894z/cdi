function [ proj_geom ] = isosurface_ASTRA_parallel3Dvec( obj, params )

% function isosurface_phase_surface_raw( obj, params )
         
                        
    amps    = abs( obj );
    phasesp = angle( obj ); %/phsf;%.*(amps>=0.15)/phsf; %p1a1 and p1a2

    if ~isfield( params, 'y_coord_length' ), params.y_coord_length = 1 : size( obj, 1 ); end
    if ~isfield( params, 'x_coord_length' ), params.x_coord_length = 1 : size( obj, 2 ); end
    if ~isfield( params, 'z_coord_length' ), params.z_coord_length = 1 : size( obj, 3 ); end
    
    if ~isfield( params, 'y_lim' ), params.y_lim = [ 1, size( obj, 1 ) ]; end
    if ~isfield( params, 'x_lim' ), params.x_lim = [ 1, size( obj, 2 ) ]; end
    if ~isfield( params, 'z_lim' ), params.z_lim = [ 1, size( obj, 3 ) ]; end
    
    [ x, y, z ] = meshgrid( params.x_coord_length, ...
                            params.y_coord_length, ...
                            params.z_coord_length );
    
    hold on;

    p = patch( isosurface( x, y, z, amps, params.isoslvl ));

    isonormals( x, y, z, amps, p );
    isocolors( x, y, z, phasesp, p ); %for phases
    set( p, 'FaceColor', 'interp', 'EdgeColor', 'none' )

    p2 = patch( isocaps( x, y, z, amps, params.isoslvl ), 'FaceColor', 'interp', 'EdgeColor','none' );
    set( p2, 'FaceColor', 'interp', 'EdgeColor', 'none' )

    ylabel( '+y', 'fontweight', 'bold', 'fontsize', 16 ); ylim( params.y_lim )
    xlabel( '+x', 'fontweight', 'bold', 'fontsize', 16 ); xlim( params.x_lim )
    zlabel( '+z', 'fontweight', 'bold', 'fontsize', 16 ); zlim( params.z_lim )
    
    view( [ 45, 30 ]);

    camlight;
    lighting gouraud;
    grid on

    alpha( params.alphalvl ) 

    hold off;

    
    
end

