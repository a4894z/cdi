function [ rs, rot ] = scanpositions_update_2DTPA_rotation_v2( rs, psi, T0, phi, spos_opt, sol_GPU, expt )

Nspos = size( rs, 1 );
Nscpm = size( phi, 3 );

% rot_search = gpuArray( linspace( -5, 5, 31 ));
rot_search = gpuArray( [ 0, 12 ]);

L = gpuArray( zeros( 1, length( rot_search ), 'single' ) );

for aa = 1 : length( rot_search )

    rotT = gpuArray( [ [ cosd( rot_search( aa ) ), -sind( rot_search( aa ) ) ]; ...
                       [ sind( rot_search( aa ) ), +cosd( rot_search( aa ) ) ] ] );     
             
    rs_rot = transpose( rotT * transpose( rs ));
            
    
    ind = get_indices_2Dframes( rs_rot, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi_rot = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;

%     ind_offset = ind + uint32( sol_GPU.samrc * ( 0 : 1 : ( sol_GPU.Nspos - 1 ) ));
%     
%     
% %          [ TF1 ] = rPIEupdate_batch_2DTPA_sample(  psi_rot,        ...
% %                                                    sol_GPU.TFvec,      ...
% %                                                    sol_GPU.phi,        ...
% %                                                    ind_offset, ...
% %                                                    sol_GPU.rc,         ...
% %                                                    sol_GPU.Nspos,      ...
% %                                                    sol_GPU.rPIE_alpha_T );  
%     
%     
%          [ TF2 ] = rPIEupdate_batch_2DTPA_sample(  psi_rot,        ...
%                                                    sol_GPU.TFvec,      ...
%                                                    sol_GPU.phi,        ...
%                                                    sol_GPU.ind_offset, ...
%                                                    sol_GPU.rc,         ...
%                                                    sol_GPU.Nspos,      ...
%                                                    sol_GPU.rPIE_alpha_T ); 
%                                                
%     
%          [ TF3 ] = rPIEupdate_batch_2DTPA_sample(  psi,        ...
%                                                    sol_GPU.TFvec,      ...
%                                                    sol_GPU.phi,        ...
%                                                    ind_offset, ...
%                                                    sol_GPU.rc,         ...
%                                                    sol_GPU.Nspos,      ...
%                                                    sol_GPU.rPIE_alpha_T ); 
%                                                
%                                                
%     TF0_mat = reshape( T0, spos_opt.szTF );
% %     TF1_mat = reshape( TF1, spos_opt.szTF );
%     TF2_mat = reshape( TF2, spos_opt.szTF );
%     TF3_mat = reshape( TF3, spos_opt.szTF );
%     
%     figure; imagescHSV(TF0_mat)
% %     figure; imagescHSV(TF1_mat)
%     figure; imagescHSV(TF2_mat)
%     figure; imagescHSV(TF3_mat)
%     
%     figure; imagesc(abs(TF0_mat))
%     figure; imagesc(abs(TF2_mat))
%     figure; imagesc(abs(TF3_mat))
    
    
    L( aa ) = sum( sum( sum( sum( abs( psi - psi_rot ) .^ 2 )))) / ( Nspos * spos_opt.rc * Nscpm );


end



[ ~, II ] = min( L );

rot = rot_search( II );

rotT = gpuArray( [ [ cosd( rot ), -sind( rot ) ] ; ...
                   [ sind( rot ), +cosd( rot ) ] ] );  
                   
rs = transpose( rotT * transpose( rs ));




figure;
plot( rot_search, L )

%{
    
    
    GPU.sample_sposview_indices = get_indices_2Dframes( GPU.batch_rs, sol.GPU.samsz, sol.GPU.vs_r, sol.GPU.vs_c ); 

    GPU.batch_N = gpuArray( length( GPU.batch_indx ) );

    GPU.Nspos = gpuArray( GPU.batch_N );
    GPU.Nscpm = gpuArray( sol.probe.scpm.N );

    %========

    GPU.ind        = uint32( GPU.sample_sposview_indices ); 
    GPU.ind_offset = uint32( GPU.samrc * ( 0 : 1 : ( GPU.Nspos - 1 ) ));
    GPU.ind_offset = GPU.ind + GPU.ind_offset;

%}


end

