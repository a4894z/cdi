function imagesc_abs_phs( X, opt )
 
if ~exist( 'opt', 'var' ), opt = struct; end

if ~isfield( opt, 'daspect' ), opt.daspect = true; end
    
ax1 = subplot( 121 );
imagesc( abs( X ))
colormap( ax1, gray )
if opt.daspect == true, daspect( [ 1, 1, 1 ] ); end 
    
ax2 = subplot( 122 );
imagesc( angle( X ))
colormap( ax2, hsv )
if opt.daspect == true, daspect( [ 1, 1, 1 ] ); end 
