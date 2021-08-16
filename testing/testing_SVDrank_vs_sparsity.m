




sol.sparse.pct         = single( 0.15 );
sol.sparse.lvl         = round( sol.sample.sz.rc * sol.sparse.pct );

sol.sparse.threshtype  = 'h';
sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 8 };
    
[ TFs ] = sparseFDxy_update_2Dsample( expt.sample.TF, sol.sparse );

figure; imagesc(abs(TFs)); daspect([1 1 1]); colormap turbo; title('15% hard thresh')


sol.sparse.threshtype  = 's';
sol.sparse.threshname  = sol.sparse.type_sparse2DFDxy{ 8 };
    
[ TFs ] = sparseFDxy_update_2Dsample( expt.sample.TF, sol.sparse );

figure; imagesc(abs(TFs)); daspect([1 1 1]); colormap turbo; title('15% soft thresh')






[ U, S, V ] = svd( expt.sample.TF );


S1 = 0 * S;
S1( :, 1 ) = S( :, 1 );

% figure; imagesc( abs( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 1st Singular Value')
figure; imagescHSV( ( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 1st Singular Value')



S1 = 0 * S;
S1( :, 2 ) = S( :, 2 );

% figure; imagesc( abs( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 2nd Singular Value')
figure; imagescHSV( ( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 2nd Singular Value')

S1 = 0 * S;
S1( :, 3 ) = S( :, 3 );

% figure; imagesc( abs( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 3rd Singular Value')
figure; imagescHSV( ( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 3rd Singular Value')

S1 = 0 * S;
S1( :, 4 ) = S( :, 4 );

% figure; imagesc( abs( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 4th Singular Value')
figure; imagescHSV( ( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 4th Singular Value')


S1 = 0 * S;
S1( :, 5 ) = S( :, 5 );

% figure; imagesc( abs( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 5th Singular Value')
figure; imagescHSV( ( U * S1 * V' )); daspect([1 1 1]); colormap turbo; title('Only 5th Singular Value')
















return




S( :, round( end * sol.sparse.pct ) : end ) = 0;

figure; imagesc( abs( U * S * V' )); daspect([1 1 1]); colormap turbo; title('top 15% singular values')














return