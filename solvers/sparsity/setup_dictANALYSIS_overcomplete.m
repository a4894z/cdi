function [ aDL ] = setup_dictANALYSIS_overcomplete( sz, expt )

%==================================================================================================

% ANALYSIS MODEL
%    Na x Np        Na x Ni                 Ni x Np
% || C_a      -     D_a                     R        ||^2
% || C_a      -     [ da_1 da_2 ... da_n ]  R        ||^2
% 
% R    is Ni x Np                            
% D_a   is Na x Ni
% C_a   is Na x Np
% 
% D_a operates on the Np image patches to create a sparse representation C_a of the image patches
% Na --> # of dictionary image feature "atoms" ??????????????????
% Ni --> # of rows for image patch vector  
% Np --> # of patches in training set
%
% aDL.wx = round( sz( 2 ) * 0.1 );          % patch size # cols
% aDL.wy = round( sz( 1 ) * 0.1 );          % patch size # rows  
%
% number of variables in training set patches ( dimension of data to be sparse coded, Np is number of patches, training set is Ni x Np ):
% aDL.Ni = aDL.wx * aDL.wy;                    
%
% number of atoms in the dictionary ( dictionary matrix D is Na x Ni ), ( image patches matrix P is Ni x Np ), ( sparse representation matrix is Na x Np )
% aDL.Na = round( 3 * aDL.Ni );     
%
% target sparsity ( number of nonzeros in each column of the sparse coding matrix )
% aDL.Cnnz = round( aDL.Na * 0.01 );              
%
% overlap for training set patches ( ??? prevent/reduce block artifacts)
% aDL.ol = 0/4;    

%==================================================================================================
%--------------------------------------- DEFAULTS -------------------------------------------------
%==================================================================================================

aDL.wx = 16 + 0 * round( sz( 2 ) * 1.0 );          
aDL.wy = 16 + 0 * round( sz( 1 ) * 1.0 );          
aDL.Ni = aDL.wx * aDL.wy;      

%==================================================================================================
%----------------------------------- ANALYSIS DICTIONARY INIT--------------------------------------
%==================================================================================================

% % RANDOM MATRIX INIT
% 
% aDL.Na = round( 3 * aDL.Ni );
% 
% aDL.D = rand( aDL.Na, aDL.Ni, 'single' ) + 1 * 1i * rand( aDL.Na, aDL.Ni, 'single' );
% 
% % % rescale WHOLE MATRIX to be unit norm:
% % aDL.D = aDL.D / norm( aDL.D, 'fro' );   
% 
% % % rescale rows to be unit norm:
% % aDL.D = aDL.D ./ repmat( sqrt( sum( abs( aDL.D ) .^ 2, 2 )), [ 1, aDL.Ni ] );   
% 
% % % rescale so that each analysis row has max value of one
% % aDL.D = aDL.D ./ repmat( max( abs( aDL.D ) , [], 2 ), [ 1, aDL.Ni ] );   
% 
% % rescale so that WHOLE MATRIX has max value of one
% aDL.D = aDL.D / max( abs( aDL.D(:) ) );   

%================================================
% USE A FORWARD DIFFERENCE MATRIX AS INIT:

e = ones( aDL.Ni, 1 );

Dax = spdiags( e, 0, aDL.Ni, aDL.Ni );
Dax = Dax + spdiags( -e, +aDL.wy, aDL.Ni, aDL.Ni );
Dax = Dax + spdiags( -e, aDL.wy - aDL.Ni, aDL.Ni, aDL.Ni );
% figure; spy(Dax)

Day = spdiags( e, 0, aDL.wy, aDL.wy );
Day = Day + spdiags( -e, +1, aDL.wy, aDL.wy );
Day = Day + spdiags( -e, 1 - aDL.wy, aDL.wy, aDL.wy );
ACell = repmat({Day}, 1, aDL.wx);
Day = blkdiag(ACell{:});
% figure; spy(Day)

Da = single( full( Dax ));
Da = [ Da; single( full( Day )) ];




% aDL.mean0 = logical( 0 ); %#ok<LOGL>
% aDL.ol = 0/3;  aDL.patch = 'raster';
% [ Ya, aDL ] = image2trainingsetpatches( expt.phiT, aDL, sz );
% 
% figure; imagesc( abs(reshape( Ya( :, 1 ), [aDL.wy, aDL.wx] )))
% Cax = Dax * double( Ya( :, 1 )); figure; imagesc( abs( reshape( Cax, [aDL.wy, aDL.wx] )))
% Cay = Day * double( Ya( :, 1 )); figure; imagesc( abs( reshape( Cay, [aDL.wy, aDL.wx] )))



% [ pinvDax ] = pinv_sparse( Dax );
% 
% Yx = pinvDax * Cax; 
% 
% figure;
% imagesc( abs( reshape( Yx, [aDL.wy, aDL.wx] )))
% 


aDL.D = full( Da );
aDL.D0 = aDL.D;
aDL.Na = size( aDL.D, 1 );

aDL.mean0 = logical( 0 ); %#ok<LOGL>
aDL.ol = 0/2;  aDL.patch = 'raster';
[ IpT, aDL ] = image2trainingsetpatches( expt.phiT, aDL, sz );

% % number of nonzeros to keep in the sparse synthesis code matrix columns
% aDL.sparse_lvl = 1.0;
% aDL.Cnnz = 0 + 1 * round( aDL.Na * aDL.sparse_lvl );              
% if aDL.Cnnz <= 0, aDL.Cnnz = 1; end
% aDL.C = aDL_sparsecode_update( IpT, aDL );


% v3 = inv( aDL.D );



aDL.C = aDL.D * IpT;  

[ v1 ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sz );
[ v2 ] = trainingsetpatches2image( aDL.C \ aDL.D, aDL, sz );

[U, S, V ] = svd( aDL.D );



figure; imagesc( abs( expt.phiT ),[0, 40]); daspect([1 1 1]); colormap jet
figure; imagesc( abs( v1 ),[0, 40]); daspect([1 1 1]); colormap jet
figure; imagesc( abs( v2 )); daspect([1 1 1]); colormap jet


figure; imagesc( abs(v1 - expt.phiT) )







aDL.D = Da;
aDL.D0 = aDL.D;

aDL.Na = size( aDL.D, 1 );

%==================================================================================================
% INITIAL ANALYSIS SPARSE REPRESENTATION

% aDL.sparse_lvl = 0.05;
% aDL.Cnnz = 0 + 1 * round( aDL.Na * aDL.sparse_lvl );              
% if aDL.Cnnz <= 0, aDL.Cnnz = 1; end
%     
% %=======================
% 
% aDL.C = rand( aDL.Na, aDL.Np, 'single' ) + 1i * rand( aDL.Na, aDL.Np, 'single' );
% % aDL.C = aDL.D * Ya;
% 
% %=======================
% 
% % determine thresholding level for each column of Ca ( size Na x Np ):
% abs_Ca = abs( aDL.C );
% temp1 = sort( abs_Ca, 'descend' );
% sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    
% 
% % hard thresholding
% aDL.C = aDL.C .* ( abs_Ca >= sel );   

%==================================================================================================
% BASED ON PARAMS DEFINED ABOVE, DETERMINE HOW MANY IMAGE PATCHES 

aDL.mean0 = logical( 0 ); %#ok<LOGL>

aDL.ol = 1/2;  aDL.patch = 'raster';
% aDL.Np = 100; aDL.patch = 'random';
% 
[ ~, aDL ] = image2trainingsetpatches( [], aDL, sz  );
% [ Ya, aDL ] = image2trainingsetpatches( imresize( expt.sample.TF, sz ), aDL, sz );

%==================================================================================================

aDL = orderfields( aDL );

%==================================================================================================
