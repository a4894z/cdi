function [ aDL ] = setup_analysis_dictionarylearning( sz )

%==================================================================================================

% ANALYSIS MODEL
%    Na x Np        Na x Ni                 Ni x Np
% || C_a      -     D_a                     Ip        ||^2
% || C_a      -     [ da_1 da_2 ... da_n ]  Ip        ||^2
% 
% Ip    is Ni x Np                            
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

aDL.wx = round( sz( 2 ) * 1.0 );          
aDL.wy = round( sz( 1 ) * 1.0 );          
aDL.Ni = aDL.wx * aDL.wy;      

aDL.Na = round( 3 * aDL.Ni );

aDL.ol = 1/3;    
aDL.patch = 'raster';
% aDL.Np = 100; aDL.patch = 'random';

aDL.mean0 = false;

%==================================================================================================
% BASED ON PARAMS DEFINED ABOVE, DETERMINE HOW MANY IMAGE PATCHES 

[ ~, aDL ] = image2trainingsetpatches( [], aDL, sz  );

%==================================================================================================
%----------------------------------- ANALYSIS DICTIONARY INIT--------------------------------------
%==================================================================================================

% % RANDOM MATRIX INIT
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
% aDL.D = ;
% aDL.D = ;
% Da * Ya;

% tmp0 = Dax * double( Ya( :, 288 ));figure; imagesc( abs( reshape( tmp0, [aDL.wy, aDL.wx] )))
% tmp0 = Day * double( Ya( :, 288 ));figure; imagesc( abs( reshape( tmp0, [aDL.wy, aDL.wx] )))
% figure; imagesc( abs(reshape( Ya( :, 288 ), [aDL.wy, aDL.wx] )))

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

aDL = orderfields( aDL );

%==================================================================================================
