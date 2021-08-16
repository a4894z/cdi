function [ sDL ] = setup_dictSYNTHESIS_overcomplete( sz, expt )

%==================================================================================================

% SYNTHESIS MODEL
%    Ni x Np         Ni x Na                Na x Np
% || R      -       D_s                    C_s       ||^2
% || R      -       [ ds_1 ds_2 ... ds_p ] C_s       ||^2  
% 
% R    is Ni x Np                            
% D_s   is Ni x Na
% C_s   is Na x Np
% 
% C_s selects out as few possible dictionary image feature "atoms" as possible
% Na --> # of dictionary image feature "atoms"
% Ni --> # of rows for image patch vector  
% Np --> # of patches in training set
% 
% sDL.wx = round( sz( 2 ) * 0.07 );          % patch size # cols
% sDL.wy = round( sz( 1 ) * 0.07 );          % patch size # rows  
% 
% number of variables in training set patches ( dimension of data to be sparse coded, Np is number of patches, training set is Ni x Np )
% sDL.Ni = sDL.wx * sDL.wy;     
% 
% number of atoms in the dictionary ( dictionary matrix D is Ni x Na ), ( sparse representation matrix C_s is Na x Np ), ( image patches matrix is Ni x Np )
% sDL.Na = round( 4 * sDL.Ni );      
% 
% target sparsity ( number of nonzeros in each column of the sparse coding matrix )
% sDL.Cnnz = round( sDL.Na * 0.02 );       
% 
% overlap for training set patches ( ??? prevent/reduce block artifacts )
% sDL.ol = 0/6;     

%==================================================================================================

sDL.is_orthog = 'orthog_no';

%==================================================================================================
% IMAGE PATCH SIZE

sDL.wx = 10 + 0 * round( sz( 2 ) * 0.025 );          
sDL.wy = 10 + 0 * round( sz( 1 ) * 0.025 );          
sDL.Ni = sDL.wx * sDL.wy;        

% INIT SYNTHESIS DICTIONARY/SPARSE CODES
sDL.D = [];
sDL.C = [];

%==================================================================================================
% USE PATCHES FROM EXISTING IMAGES AS DICTIONARY "WORDS", DICTIONARY IS Ni x Na

sDL.mean0 = logical( 0 ); %#ok<LOGL>                    % subtract off the mean of image feature patches?
% sDL.D_unitnorm = logical( 1 ); %#ok<*LOGL>
% sDL.D_rescale_maxval1 = logical( 0 );
% sDL.D_project_maxval1 = logical( 0 );

% sDL.Np = sDL.Ni; sDL.patch = 'random';                  % dictionary will be complete, not over complete
% sDL.Np = 300; sDL.patch = 'random';                     % dictionary can be whatever ( complete, over/under )
sDL.ol = 1/2; sDL.patch = 'raster';                   

img_dir = expt.paths.rimgdata;

% use patches from images to define columns of dictionary atoms, say 0.2 from each, using random locations 

isDL_im = [];
isDL_im{ end + 1 } = 'sample/tulips.ppm';
isDL_im{ end + 1 } = 'sample/btvn.ppm';
isDL_im{ end + 1 } = 'sample/cells_string_aS.ppm';
isDL_im{ end + 1 } = 'sample/fractal1c.ppm';
isDL_im{ end + 1 } = 'sample/leaf2.ppm';
isDL_im{ end + 1 } = 'sample/spiral2.ppm';
isDL_im{ end + 1 } = 'sample/xy.ppm';
isDL_im{ end + 1 } = 'sample/fractal1.ppm';
isDL_im{ end + 1 } = 'sample/prozorov.ppm';
isDL_im{ end + 1 } = 'sample/xy3.ppm';
isDL_im{ end + 1 } = 'sample/baboon.ppm';
isDL_im{ end + 1 } = 'sample/cells.ppm';
isDL_im{ end + 1 } = 'sample/hsv.ppm';
isDL_im{ end + 1 } = 'sample/ising2.ppm';
isDL_im{ end + 1 } = 'sample/phase_spiral.ppm';
isDL_im{ end + 1 } = 'sample/SheppLogan2.ppm';
% isDL_im{ end + 1 } = 'sample/siemens-spokeB.ppm';
isDL_im{ end + 1 } = 'sample/sw_a.ppm';
isDL_im{ end + 1 } = 'sample/sw_b.ppm';

Ds = [];

for ii = 1 : length( isDL_im )
    
%     d = imread( [ img_dir, isDL_im{ ii }]);
%     d = rgb2hsv( d );
%     d = d( :, :, 3 ) .* exp( 1i * 1 * pi * d( :, :, 1 ));      % brightness --> abs, hue --> phase
    
    [ d ] = image2complex( [ img_dir, isDL_im{ ii }], 3, 1 );
    
    d = single( imresize( d, [ sz( 1 ), sz( 2 ) ] ));
        
    d = modulus_limits_scale( d, [ 0, 1 ] );
    d = phase_limits_scale( d, [ -pi, +pi ] );
    
    [ D( :, :, ii ), DL{ ii } ] = image2trainingsetpatches( d, sDL, sz );
    
    
%     % orthogonalize the columns:
%     [ u, s, v ] = svd( D( :, :, ii ) );
%     D( :, :, ii ) = u * v';
    
    Ds = [ Ds, D( :, :, ii ) ];
    
    
end

% ASSIGN THESE IMAGE FEATURES TO THE SYNTHESIS DICTIONARY
sDL.D = Ds;
sDL.Na = size( sDL.D, 2 );

% % RANDOMLY ASSIGN THESE IMAGE FEATURES TO THE SYNTHESIS DICTIONARY
% sDL.Na = 1 * sDL.Ni;
% sDL.D = Ds( :, randperm( size( Ds, 2 ), sDL.Na ));

% display the initial dictionary.
figure; plot_dictionnary( sDL, [ sDL.wy, sDL.wx ] );
close all;










% tilde_S = sDL.D;
% [ U, ~ ] = eig( tilde_S' * tilde_S, 'nobalance', 'chol' );  % 'chol', 'qz'
% S = tilde_S * U; 
% 
% 
% Nr = 512; 
% Nc = 4; 
% A = rand( Nr, Nc ) + 1i * rand( Nr, Nc ); 
% [ U, ~ ] = eig( A' * A, 'nobalance' ); 
% S = A * U;







%==================================================================================================
%scaling parameters for synthesis dictionary 

% sDL.D = sDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( sDL.D ) .^ 2 )), [ sDL.Ni, 1 ] ));      % rescale columns ( synthesis dictionary terms ) to be unit norm:

% sDL.D = sDL.D ./ ( 1e-7 + repmat( max( abs( sDL.D )), [ sDL.Ni, 1 ] ));                   % rescale so that each synthesis dictionary component has max value of one

% tmp0 = abs( sDL.D ) > 1; sDL.D( tmp0 ) = 1 * exp( 1i * angle( sDL.D( tmp0 )));            % project so that each abs() synthesis dictionary component has max value of one

%==================================================================================================
% DEFINE HOW EXTRACT IMAGE PATCHES, IMAGE PATCHES ARE Ni x Np ( THIS DEFINES Np )

sDL.ol = 0/3; sDL.patch = 'raster';
% sDL.Np = 1000; sDL.patch = 'random';

[ ~, sDL ] = image2trainingsetpatches( [], sDL, sz  );

%==================================================================================================
% INITIAL SYNTHESIS SPARSE CODE REPRESENTATION, SPARSE CODE IS Na x Np
% 
% % sDL.C = eye( sDL.Na, sDL.Np, 'single' );
% sDL.C = rand( sDL.Na, sDL.Np, 'single' ) + 1i * rand( sDL.Na, sDL.Np, 'single' );
% 
% %=======================
% 
% abs_C = abs( sDL.C );
% thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * 0.05 );
% sDL.C = sDL.C .* ( abs_C >= thresh );  

%=======================

% sDL.sparse_lvl = 0.1;
% sDL.Cnnz = 0 + 1 * round( sDL.Na * sDL.sparse_lvl );
% if sDL.Cnnz <= 0, sDL.Cnnz = 1; end
% 
% % determine thresholding level for each column of Cs:
% abs_Cs = abs( sDL.C );
% temp1 = sort( abs_Cs, 'descend' );
% sel = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );  
% 
% % hard thresholding
% sDL.C = sDL.C .* ( abs_Cs >= sel );   
% 
% sDL.Cabsminmax = [ 0, 10 ];
% sDL.Cphsminmax = 1.0 * [ -pi, pi ];
% [ sDL.C ] = scale_matrixcolumns_magnitude( sDL.C, sDL.Cabsminmax );
% [ sDL.C ] = scale_matrixcolumns_phase( sDL.C, sDL.Cphsminmax );

%==================================================================================================

sDL = orderfields( sDL );

%==================================================================================================



