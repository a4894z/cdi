function [ proj_x, proj_y, proj_z ] = imagesc_3Dxyz_projections( V3D, kwargs_proj, kwargs_colormap )

if ~exist( 'kwargs_proj', 'var' ), kwargs_proj = 'sum_proj'; end
if isempty( kwargs_proj ),         kwargs_proj = 'sum_proj'; end

if ~exist( 'kwargs_colormap', 'var' ), kwargs_colormap = bone; end
if isempty( kwargs_colormap ),         kwargs_colormap = bone; end

%========

switch kwargs_proj
    
    case 'sum_proj'
        
        proj_x = squeeze( sum( V3D, 2 ) / numel( size( V3D, 2 )));
        proj_y = squeeze( sum( V3D, 1 ) / numel( size( V3D, 1 )));
        proj_z = squeeze( sum( V3D, 3 ) / numel( size( V3D, 3 )));
        
        a1 = subplot( 1, 3, 1 );
        imagesc( proj_y )
        daspect( [ 1, 1, 1 ])
        colormap( a1, kwargs_colormap );
        colorbar
        title( 'proj\_y' )
      
        a2 = subplot( 1, 3, 2 );
        imagesc( proj_x )
        daspect( [ 1, 1, 1 ])
        colormap( a2, kwargs_colormap );
        colorbar
        title( 'proj\_x' )
        
        a3 = subplot( 1, 3, 3 );
        imagesc( proj_z )
        daspect( [ 1, 1, 1 ])
        colormap( a3, kwargs_colormap );
        colorbar
        title( 'proj\_z' )
                
    case 'prod_proj'
        
        proj_x = prod( V3D, 2 );
        proj_y = prod( V3D, 1 );
        proj_z = prod( V3D, 3 );

        a1 = subplot( 1, 3, 1 );
        imagesc( proj_y )
        daspect( [ 1, 1, 1 ])
        colormap( a1, kwargs_colormap );
        colorbar
        title( 'proj\_y' )
        
        a2 = subplot( 1, 3, 2 );
        imagesc( proj_x )
        daspect( [ 1, 1, 1 ])
        colormap( a2, kwargs_colormap );
        colorbar
        title( 'proj\_x' )
        
        a3 = subplot( 1, 3, 3 );
        imagesc( proj_z )
        daspect( [ 1, 1, 1 ])
        colormap( a3, kwargs_colormap );
        colorbar
        title( 'proj\_z' )
        
    otherwise
        
        error( "Incorrect 'kwargs_proj'; use 'sum_proj' or 'prod_proj' " )
        
end



end

