function imagesc_3Dxyz_montage( V3D, slice_montage )

kk = 1;

%========

switch slice_montage.rcp
    
    case 'rows'
        
        slice_montage.skip         = floor( size( V3D, 1 ) / prod( slice_montage.nrpw_ncols ));
        slice_montage.range        = 1 : slice_montage.skip : size( V3D, 1 );

        while length( slice_montage.range ) > prod( slice_montage.nrpw_ncols )
            
            slice_montage.range( end ) = [];
            
        end

        for slc = slice_montage.range
        
            a1 = subplot( slice_montage.nrpw_ncols( 1 ), slice_montage.nrpw_ncols( 2 ), kk );
            imagesc( squeeze( V3D( slc, :, : )), slice_montage.lims )
            daspect( [ 1, 1, 1 ])
            colormap( a1, slice_montage.colormap );
            colorbar
            
            kk = kk + 1;
        
        end

    case 'cols'
        
        slice_montage.skip         = floor( size( V3D, 2 ) / prod( slice_montage.nrpw_ncols ));
        slice_montage.range        = 1 : slice_montage.skip : size( V3D, 2 );
        
        while length( slice_montage.range ) > prod( slice_montage.nrpw_ncols )
            
            slice_montage.range( end ) = [];
            
        end

        for slc = slice_montage.range
        
            a1 = subplot( slice_montage.nrpw_ncols( 1 ), slice_montage.nrpw_ncols( 2 ), kk );
            imagesc( squeeze( V3D( :, slc, : )), slice_montage.lims )
            daspect( [ 1, 1, 1 ])
            colormap( a1, slice_montage.colormap );
            colorbar
            
            kk = kk + 1;
        
        end
        
    case 'pages'
        
        slice_montage.skip         = floor( size( V3D, 3 ) / prod( slice_montage.nrpw_ncols ));
        slice_montage.range        = 1 : slice_montage.skip : size( V3D, 3 );

        while length( slice_montage.range ) > prod( slice_montage.nrpw_ncols )
            
            slice_montage.range( end ) = [];
            
        end

        for slc = slice_montage.range
        
            a1 = subplot( slice_montage.nrpw_ncols( 1 ), slice_montage.nrpw_ncols( 2 ), kk );
            imagesc( squeeze( V3D( :, :, slc )), slice_montage.lims )
            daspect( [ 1, 1, 1 ])
            colormap( a1, slice_montage.colormap );
            colorbar
            
            kk = kk + 1;
        
        end
        
    otherwise
        
        error( "Incorrect 'kwargs_proj'; use 'sum_proj' or 'prod_proj' " )
        
end

%========

% switch slice_montage.rcp
%     
%     case 'rows'
%         
%         title( 'row slices' )
% 
%     case 'cols'
%         
%         title( 'column slices' )
%         
%     case 'pages'
%         
%         title( 'page slices' )   
%         
%     otherwise
%         
%         error( "Incorrect 'slice_montage.rcp'; use 'rows' or 'cols' or 'pages' " )
%         
% end


