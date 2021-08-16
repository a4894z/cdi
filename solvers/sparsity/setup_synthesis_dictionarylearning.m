function [ sDL ] = setup_synthesis_dictionarylearning( sz, expt )

%==================================================================================================

% SYNTHESIS MODEL
%    Ni x Np         Ni x Na                Na x Np
% || Ip      -       D_s                    C_s       ||^2
% || Ip      -       [ ds_1 ds_2 ... ds_p ] C_s       ||^2  
% 
% Ip    is Ni x Np                            
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
% IMAGE PATCH SIZE

sDL.wx = 0 + 0 * round( sz( 2 ) * 1.00 );          
sDL.wy = 0 + 0 * round( sz( 1 ) * 1.00 );          
sDL.Ni = sDL.wx * sDL.wy;        

% INIT SYNTHESIS DICTIONARY/SPARSE CODES
sDL.D = [];
sDL.C = [];

%==================================================================================================
% SYNTHESIS DICTIONARY INIT USING RANDOM NUMBERS

%{

sDL.Na = round( 1 * sDL.Ni );
sDL.D = rand( sDL.Ni, sDL.Na, 'single' ) + 1 * 1i * rand( sDL.Ni, sDL.Na, 'single' );

%}

%==================================================================================================
% USE PATCHES FROM EXISTING IMAGES AS DICTIONARY "WORDS", DICTIONARY IS Ni x Na

sDL.mean0 = logical( 1 ); %#ok<LOGL>                    % subtract off the mean of image feature patches?
% sDL.D_unitnorm = logical( 1 ); %#ok<*LOGL>
% sDL.D_rescale_maxval1 = logical( 0 );
% sDL.D_project_maxval1 = logical( 0 );

% sDL.Np = sDL.Ni; sDL.patch = 'random';                  % dictionary will be complete, not over complete
sDL.Np = 100; sDL.patch = 'random';                     % dictionary can be whatever ( complete, over/under )
% sDL.ol = 0/3; sDL.patch = 'raster';                   

img_dir = '~/Documents/MATLAB/Code/data/simulated/sample/';

% use patches from images to define columns of dictionary atoms, say 0.2 from each, using random locations 

isDL_im = [];
isDL_im{ end + 1 } = 'tulips.ppm';
% isDL_im{ end + 1 } = 'btvn.ppm';
isDL_im{ end + 1 } = 'cells_string_aS.ppm';
isDL_im{ end + 1 } = 'fractal1c.ppm';
isDL_im{ end + 1 } = 'leaf2.ppm';
isDL_im{ end + 1 } = 'spiral2.ppm';
isDL_im{ end + 1 } = 'xy.ppm';
isDL_im{ end + 1 } = 'fractal1.ppm';
isDL_im{ end + 1 } = 'prozorov.ppm';
isDL_im{ end + 1 } = 'xy3.ppm';
isDL_im{ end + 1 } = 'baboon.ppm';
isDL_im{ end + 1 } = 'cells.ppm';
isDL_im{ end + 1 } = 'hsv.ppm';
isDL_im{ end + 1 } = 'ising2.ppm';
isDL_im{ end + 1 } = 'phase_spiral.ppm';
isDL_im{ end + 1 } = 'SheppLogan2.ppm';
isDL_im{ end + 1 } = 'siemens-spokeB.ppm';
isDL_im{ end + 1 } = 'sw_a.ppm';
isDL_im{ end + 1 } = 'sw_b.ppm';

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
    Ds = [ Ds, D( :, :, ii ) ];
    
end

% % ASSIGN THESE IMAGE FEATURES TO THE SYNTHESIS DICTIONARY
% sDL.D = Ds;
% sDL.Na = size( sDL.D, 2 );

% RANDOMLY ASSIGN THESE IMAGE FEATURES TO THE SYNTHESIS DICTIONARY
sDL.Na = 1 * sDL.Ni;
sDL.D = Ds( :, randperm( size( Ds, 2 ), sDL.Na ));

% % display the initial dictionary.
% figure; plot_dictionnary( sDL, [ 16, 16 ] );
% close all;

%==================================================================================================

%{

% IP.Np = sDL.Ni; IP.patch = 'random';
% IP.Np = 200; IP.patch = 'random';
IP.ol = 0/3; IP.patch = 'raster';
IP.Ni = sDL.Ni;
IP.wx = sDL.wx;
IP.wy = sDL.wy;

IP.mean0 = 0;

[ D, ~ ] = image2trainingsetpatches( expt.sample.TF, IP, expt.sample.sz.sz );

% ASSIGN THESE IMAGE FEATURES TO THE SYNTHESIS DICTIONARY
sDL.D = [ sDL.D, D ];
sDL.Na = size( sDL.D, 2 );

% figure; plot_dictionnary( sDL, [ 8, 8 ] );
% close all;
% 5;
% 
% figure; imagesc( abs( sDL.D' * sDL.D))

%}

%==================================================================================================
%scaling parameters for synthesis dictionary 

% figure; plot_dictionnary( sDL, [ 16, 16 ] );

% sDL.D = sDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( sDL.D ) .^ 2 )), [ sDL.Ni, 1 ] ));      % rescale columns ( synthesis dictionary terms ) to be unit norm:

% sDL.D = sDL.D ./ ( 1e-7 + repmat( max( abs( sDL.D )), [ sDL.Ni, 1 ] ));                   % rescale so that each synthesis dictionary component has max value of one

% tmp0 = abs( sDL.D ) > 1; sDL.D( tmp0 ) = 1 * exp( 1i * angle( sDL.D( tmp0 )));            % project so that each abs() synthesis dictionary component has max value of one

% sDL.D = sDL.D' * sDL.D;                                                                   % make Hermitian

[ u, s, v ] = svd( sDL.D );
% s = eye( sDL.Ni, sDL.Na );
% sDL.D = u * s * v';
sDL.D = u * v';                                                                           % orthogonalize the columns


% sDL.Dabsminmax = [ 0, 1 ];
% sDL.Dphsminmax = 0.5 * [ -pi, pi ];
% [ sDL.D ] = scale_matrixcolumns_phase( sDL.D, sDL.Dphsminmax );
% [ sDL.D ] = scale_matrixcolumns_magnitude( sDL.D, sDL.Dabsminmax );

% % figure; imagesc( abs( sDL.D ))
% % figure; imagesc( abs( sDL.D' * sDL.D ))
% figure; plot_dictionnary( sDL, [ 16, 16 ] );
% close all;

%==================================================================================================
% DEFINE HOW TO EXTRACT IMAGE PATCHES, IMAGE PATCHES ARE Ni x Np ( THIS DEFINES Np )

sDL.ol = 1/2; sDL.patch = 'raster';
% sDL.Np = 1000; sDL.patch = 'random';

[ ~, sDL ] = image2trainingsetpatches( [], sDL, sz  );

%==================================================================================================
% INITIAL SYNTHESIS SPARSE CODE REPRESENTATION, SPARSE CODE IS Na x Np

% sDL.C = eye( sDL.Na, sDL.Np, 'single' );
sDL.C = rand( sDL.Na, sDL.Np, 'single' ) + 1i * rand( sDL.Na, sDL.Np, 'single' );

%=======================

abs_C = abs( sDL.C );
thresh = find_thresh_from_sparsitylevel( abs_C, numel( abs_C ) * 0.05 );
sDL.C = sDL.C .* ( abs_C >= thresh );  

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




%{



% DISCRETE COSINE TRANSFORM BASIS FUNCTIONS

clear; close all;

Nr = 8; Nc = 8;

[ pp, qq ] = meshgrid( 1 : Nc, 1 :  Nr );


nr_max = 8; 
nc_max = 8;

ii = 1;
jj = 1;

tmp0 = zeros( Nr, Nc, nr_max, nc_max, 'single' );
tmp1 = zeros( Nr * Nc, nr_max * nc_max, 'single' );
tmp2 = zeros( Nr * nr_max, Nc * nc_max, 'single' );

for nr = 0 : nr_max - 1
    for nc = 0 : nc_max - 1

        if nr == 0, alpha_nr = 1 / sqrt( Nr );
        else, alpha_nr = sqrt( 2 /  Nr ); end

        if nc == 0, alpha_nc = 1 / sqrt( Nc );
        else, alpha_nc = sqrt( 2 /  Nc ); end


        % B_pq = alpha_nr * alpha_nc * cos( pi * ( 2 * pp + 1 ) * nc / ( 2 * Nc )) .* cos( pi * ( 2 * qq + 1 ) * nr / ( 2 * Nr ));
        B_pq = alpha_nr * alpha_nc * cos( pi * ( 2 * (1 : Nr).' + 1 ) * nr / ( 2 * Nr )) * cos( pi * ( 2 * (1 : Nc) + 1 ) * nc / ( 2 * Nc ));
%         B_pq = alpha_nr * alpha_nc * cos( pi * ( 2 * nr + 1 ) * (1 : Nr).' / ( 2 * Nr )) * cos( pi * ( 2 * nc + 1 ) * (1 : Nc) / ( 2 * Nc ));
        
        tmp0( :, :, ii, jj ) = B_pq;
        tmp1( :, ( ii - 1 ) * nc_max + jj ) = B_pq( : );
        
        tmp2( nr * Nr + ( 1 : Nr ), nc * Nc + ( 1 : Nc ) ) = B_pq;
        
        jj = jj + 1;
        
        figure; imagesc( B_pq ); colormap gray; colorbar
        close all;
        
    end
    
    ii = ii + 1;
    jj = 1;
end

% figure; imagesc( tmp1 ); colormap gray; colorbar
% figure; imagesc( tmp2 ); colormap gray; colorbar



figure; imagesc( tmp1 ); colormap gray; colorbar
figure; imagesc( tmp1' * tmp1 ); colormap gray; colorbar

% orthogonalize the columns:
[ u, s, v ] = svd( tmp1 );
% s = eye( sDL.Ni, sDL.Na );
% sDL.D = u * s * v';
tmp1 = u * v';

figure; imagesc( tmp1' * tmp1 ); colormap gray; colorbar
figure; imagesc( tmp1 ); colormap gray; colorbar


% tmp2( 0 * Nr + ( 1 : Nr ), 0 * Nc + ( 1 : Nc ) ) = B_pq;
% tmp2( 0 * Nr + ( 1 : Nr ), 1 * Nc + ( 1 : Nc ) ) = B_pq;
% tmp2( 0 * Nr + ( 1 : Nr ), 2 * Nc + ( 1 : Nc ) ) = B_pq;
% 
% tmp2( 1 * Nr + ( 1 : Nr ), 0 * Nc + ( 1 : Nc ) ) = B_pq;
% tmp2( 1 * Nr + ( 1 : Nr ), 1 * Nc + ( 1 : Nc ) ) = B_pq;
% tmp2( 1 * Nr + ( 1 : Nr ), 2 * Nc + ( 1 : Nc ) ) = B_pq;

return


%}

