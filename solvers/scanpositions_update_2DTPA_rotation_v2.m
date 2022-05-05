function [ rs, rot ] = scanpositions_update_2DTPA_rotation_v2( rs, psi, T0, phi, measD, nDeq0, spos_opt, sol_GPU, expt )

nDeq0 = not( nDeq0 );
Nspos = size( rs, 1 );
Nscpm = size( phi, 3 );

% rot_search = gpuArray( linspace( -5, 5, 31 ));
rot_search = gpuArray( [ 12, 0, -12 ]);

L = gpuArray( zeros( 1, length( rot_search ), 'single' ) );

for aa = 1 : length( rot_search )

    rotT = gpuArray( [ [ cosd( rot_search( aa ) ), -sind( rot_search( aa ) ) ]; ...
                       [ sind( rot_search( aa ) ), +cosd( rot_search( aa ) ) ] ] );     
             
    rs_rot = transpose( rotT * transpose( rs ));
            
    
    ind = get_indices_2Dframes( rs_rot, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
    psi_rot = reshape( T0( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;

    ind_offset = ind + uint32( sol_GPU.samrc * ( 0 : 1 : ( sol_GPU.Nspos - 1 ) ));
    
    
    
    
    TF0_mat = reshape( T0, spos_opt.szTF );

    
    
%          [ TF1 ] = rPIEupdate_batch_2DTPA_sample(  psi_rot,        ...
%                                                    sol_GPU.TFvec,      ...
%                                                    sol_GPU.phi,        ...
%                                                    ind_offset, ...
%                                                    sol_GPU.rc,         ...
%                                                    sol_GPU.Nspos,      ...
%                                                    sol_GPU.rPIE_alpha_T );  
%     
%         TF1_mat = reshape( TF1, spos_opt.szTF );
    
    
         [ TF2 ] = rPIEupdate_batch_2DTPA_sample(  psi_rot,        ...
                                                   sol_GPU.TFvec,      ...
                                                   sol_GPU.phi,        ...
                                                   sol_GPU.ind_offset, ...
                                                   sol_GPU.rc,         ...
                                                   sol_GPU.Nspos,      ...
                                                   sol_GPU.rPIE_alpha_T ); 
                                               
      TF2_mat = reshape( TF2, spos_opt.szTF );                                               
    
      
%      ind = get_indices_2Dframes( rs_rot, spos_opt.szTF, spos_opt.vs_r, spos_opt.vs_c ); 
     psi_rot = reshape( TF2( ind ), [ spos_opt.sz, 1, Nspos ]) .* phi;
     psi_rot = fft( fft( fftshift( fftshift( psi_rot, 1 ), 2 ), [], 1 ), [], 2 ) / spos_opt.sqrt_rc;
      
     
     
TF0rot_a = imrotate( TF0_mat, -12 , 'nearest', 'crop' );

figure; imagescHSV(TF0_mat); title('TF0\_mat')
figure; imagescHSV(TF0rot_a); title('TF0rot\_a')
 figure; imagescHSV( expt.sample.T)




%          [ TF3 ] = rPIEupdate_batch_2DTPA_sample(  psi,        ...
%                                                    sol_GPU.TFvec,      ...
%                                                    sol_GPU.phi,        ...
%                                                    ind_offset, ...
%                                                    sol_GPU.rc,         ...
%                                                    sol_GPU.Nspos,      ...
%                                                    sol_GPU.rPIE_alpha_T ); 
%                                                
%     TF3_mat = reshape( TF3, spos_opt.szTF );
%     




    figure; imagescHSV(TF0_mat); title('TF0\_mat')
%     figure; imagescHSV(TF1_mat); title('TF1\_mat')
    figure; imagescHSV(TF2_mat); title('TF2\_mat')
%     figure; imagescHSV(TF3_mat); title('TF3\_mat')
    figure; imagescHSV( expt.sample.T)
    
%     figure; imagesc(abs(TF0_mat))
%     figure; imagesc(abs(TF2_mat))
%     figure; imagesc(abs(TF3_mat))
    
    sqrt_I_e      = sqrt( sum( abs( psi_rot ) .^ 2, 3 ));
%     I_e      = squeeze( sum( abs( psi ) .^ 2, 3 ) );

    L( aa ) = sum( sum( sum( abs( measD - nDeq0 .* sqrt_I_e  ) .^ 2, 1 ), 2 )) / ( Nspos * spos_opt.rc );
    
%     L( aa ) = sum( sum( sum( sum( abs( psi - psi_rot ) .^ 2 )))) / ( Nspos * spos_opt.rc * Nscpm );


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

