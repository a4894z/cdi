function [ v, sDL ] = proj_sparse_synthesis_learning( v, sDL, sz )

%==================================================================================================
%-------------------------------------- Extract image patches -------------------------------------
%==================================================================================================

[ Ip, sDL ] = image2trainingsetpatches( v, sDL, sz  );

%==================================================================================================
%------------------ Define initial sparse codes and dictionary if not defined ---------------------
%==================================================================================================

% if isempty( sDL.C ) || ( size( sDL.C, 1 ) ~= sDL.Na ) || ( size( sDL.C, 2 ) ~= sDL.Np )
% 
%     % sparse SYNTHESIS code init
%     sDL.C = rand( sDL.Na, sDL.Np, 'single' ) + 1i * rand( sDL.Na, sDL.Np, 'single' );
%     
%     % enforce sparsity level on sparse code vectors:
%     abs_Cs = abs( sDL.C );
%     temp1 = sort( abs_Cs, 'descend' );
%     sel = repmat( temp1( sDL.Cnnz, : ), [ size( temp1, 1 ), 1 ] );    
%     % hard thresholding
%     sDL.C = sDL.C .* ( abs_Cs >= sel );  
%    
% end
% 
% if isempty( sDL.D ) || ( size( sDL.D, 1 ) ~= sDL.Ni ) || ( size( sDL.D, 2 ) ~= sDL.Na )
%     
%     % SYNTHESIS dictionary init:
%     sDL.D = rand( sDL.Ni, sDL.Na, 'single' ) + 1i * rand( sDL.Ni, sDL.Na, 'single' );
%     
%     % columns of dictionary (atoms) have unit norm
%     sDL.D = sDL.D ./ ( 1e-7 + repmat( sqrt( sum( abs( sDL.D ) .^ 2, 1 )), [ sDL.Ni, 1 ] ));     
% 
% end

%==================================================================================================
%-------------------- Run some updates on the sparse codes and dictionary matrix ------------------
%==================================================================================================

for jj = 1 : sDL.ntot
    
%     if strcmp( sDL.quiet, 'no' )
%         fprintf( [num2str( [ jj, sDL.ntot ], '%d/%d' ),'\n'] ); end
    
    make_iteration_dots( 'Testing, Synthesis Sparse Learning ', jj, jj, 25 ); 
    
    tic;
        
    %============================
    % --- coefficient update ----
    %============================

    [ sDL.C, sDL.stepCs( jj ) ] = sparse_synthesis_learn_coeff_update( sDL.C, sDL.D, Ip, sDL );

    sDL.EC( end + 1 ) = sum( sum( abs( Ip - sDL.D * sDL.C ).^2 ));
    
    %============================
    % ---- dictionary update ----
    %============================
    
    if ( mod( jj, 1 ) == 0 ) && ( jj > 0 ) || ( jj == 2e99 )
            
        [ sDL.D, sDL.stepDs( jj ) ] = sparse_synthesis_learn_dict_update( sDL.C, sDL.D, Ip, sDL );
     
    end
    
    sDL.ED( end + 1 ) = sum( sum( abs( Ip - sDL.D * sDL.C ).^2 ));
            
    sDL.t( jj ) = toc;

end

%==================================================================================================
%----------- Using the sparse codes Cs, and dictionary components Ds, restore the image -----------
%==================================================================================================

[ v ] = trainingsetpatches2image( sDL.D * sDL.C, sDL, sz );

