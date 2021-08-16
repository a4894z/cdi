function [ v, aDL ] = proj_sparse_analysis_learning( v, aDL, sz )

%==================================================================================================
%-------------------------------------- Extract image patches -------------------------------------
%==================================================================================================

[ Ip, aDL ] = image2trainingsetpatches( v, aDL, sz );

%==================================================================================================
%------------------ Define initial sparse codes and dictionary if not defined ---------------------
%==================================================================================================

% if isempty( aDL.C ) || ( size( aDL.C, 1 ) ~= aDL.Na ) || ( size( aDL.C, 2 ) ~= aDL.m )
%     
%     % sparse ANALYSIS code init
%     aDL.C = rand( aDL.Na, aDL.m, 'single' ) + 1 * 1i * rand( aDL.Na, aDL.m, 'single' );
% 
%     % enforce sparsity level on sparse code vectors:
%     abs_Ca = abs( aDL.C );
%     temp1 = sort( abs_Ca, 'descend' );
%     sel = repmat( temp1( aDL.k, : ), [ size( temp1, 1 ), 1 ] );    
%     % hard thresholding
%     aDL.C = aDL.C .* ( abs_Ca >= sel );  
% 
% end
% 
% if isempty( aDL.D ) || ( size( aDL.D, 1 ) ~= aDL.Na ) || ( size( aDL.D, 2 ) ~= aDL.Ni )
%     
%     % ANALYSIS dictionary init
%     aDL.D = rand( aDL.Na, aDL.Ni, 'single' ) + 1 * 1i * rand( aDL.Na, aDL.Ni, 'single' );
%     
%     % rows of dictionary (atoms) have unit norm
%     aDL.D = aDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( aDL.D ) .^ 2, 2 )), [ 1, aDL.Ni ] ));     
%     
% %     % cols of dictionary have unit norm
% %    aDL.D = aDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( aDL.D ) .^ 2, 1 )), [ aDL.Na, 1 ] )); 
%   
% end

%==================================================================================================
%-------------------- Run some updates on the sparse codes and dictionary matrix ------------------
%==================================================================================================

for jj = 1 : aDL.ntot
    
%     if strcmp( aDL.quiet, 'no' )
%         fprintf( [num2str( [ jj, aDL.ntot ], '%d/%d' ),'\n'] ); end

    make_iteration_dots( 'Testing, Analysis Sparse Learning ', jj, jj, 25 ); 
    
    tic;
    
    %============================
    % --- coefficient update ----
    %============================

    [ aDL.C, ~ ] = sparse_analysis_learn_coeff_update( aDL.C, aDL.D, Ip, aDL );

    
    
    
    
    
    
    
    
    
    
    aDL.C = Da * Ip;      % minimizer for min_Ca || Da * Ip - Ca ||^2_F
   
    %==============================================================================================
    
    % determine thresholding level for each column of X:
    abs_Ca = abs( aDL.C );
    temp1 = sort( abs_Ca, 'descend' );
    sel = repmat( temp1( aDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    

    %============================================
    
%     % hard thresholding
%     Ca = Ca .* ( abs_Ca >= sel );  
    
    %============================================
    
    % soft thresholding
    aDL.C = ( abs_Ca -  1.0 * sel ) .* ( abs_Ca >= sel ) .* ( aDL.C ./ ( 1e-7 + abs_Ca ));  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %aDL.E0( end + 1 ) = sum(sum( abs( aDL.C - aDL.D * Ip ).^2 ));
    
    %{
    
    [ v ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sz );
    figure; imagesc( abs(v) )
    
    %}
    
    %============================
    % ---- dictionary update ----
    %============================
    
    if ( mod( jj, 1 ) == 0 )
        
        [ aDL.D, aDL.stepDa( jj ) ] = sparse_analysis_learn_dict_update( aDL.C, aDL.D, Ip, aDL );

    end
    
    aDL.E0( end + 1 ) = sum(sum( abs( aDL.C - aDL.D * Ip ).^2 ));
    
%     [ vi ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sz );
%     aDL.E0( end + 1 ) = sum(sum( abs( vi - v ).^2 ));
    
    aDL.t( jj ) = toc;

end

%==================================================================================================
%----------- Using the sparse codes Ca, and dictionary components Da, restore the image -----------
%==================================================================================================

[ v ] = trainingsetpatches2image( pinv( aDL.D ) * aDL.C, aDL, sz );

