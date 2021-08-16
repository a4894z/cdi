function [ phi ] = RAARupdate_exwv_sparseFDxy_meas( phi, sFDxy, sol, expt )

Pi_s_phi = phi;
Pi_m_R_s_phi = phi;
Pi_m_phi = phi;

for ss = sol.spos.updateorder
    
    rs = +sol.spos.rs( ss, : );             % !!!!!!!!!!!
%     rs = -sol.spos.rs( ss, : );
%     rs( 1 ) = -rs( 1 );

    % 0.5 * bbeta * ( R_m[ R_s[ phi ]] + phi ) + ( 1 - bbeta ) * Pi_m[ phi ]

    %==============================================================================================

    for pp = 1 : sol.probe.scpm.N 
    

        [ phiFD ] = edgedetect_FDxy( phi( :, :, pp, ss ), sFDxy );

        abs_phiFD.y = abs( phiFD.x );
        abs_phiFD.x = abs( phiFD.y );

        %=======================

        if strcmp( sFDxy.nzmask, 'isotropic' )

            tmp0 = sqrt( abs_phiFD.y .^ 2 + abs_phiFD.x .^ 2 );
            thresh = find_thresh_from_sparsitylevel( tmp0, sFDxy.sparselvl );

            if strcmp( sFDxy.thresh, 'hard' )

                [ phiFD.y, Wy ] = hard_shrinkage( phiFD.y, tmp0, thresh );
                [ phiFD.x, Wx ] = hard_shrinkage( phiFD.x, tmp0, thresh );

            elseif strcmp( sFDxy.thresh, 'soft' )

                [ phiFD.y, Wy ] = soft_shrinkage( phiFD.y, tmp0, thresh );
                [ phiFD.x, Wx ] = soft_shrinkage( phiFD.x, tmp0, thresh );

            else

                error( 'SPECIFY CORRECT TYPE OF THRESHOLDING (''hard'' or ''soft'' for sFDxy.thresh), EXITING )' )

            end

        elseif strcmp( sFDxy.nzmask, 'anisotropic' )

            threshr = find_thresh_from_sparsitylevel( abs_phiFD.y, sFDxy.sparselvl );
            threshc = find_thresh_from_sparsitylevel( abs_phiFD.x, sFDxy.sparselvl );

            if strcmp( sFDxy.thresh, 'hard' )

                [ phiFD.y, Wy ] = hard_shrinkage( phiFD.y, abs_phiFD.y, threshr );
                [ phiFD.x, Wx ] = hard_shrinkage( phiFD.x, abs_phiFD.x, threshc );

            elseif strcmp( sFDxy.thresh, 'soft' )

                [ phiFD.y, Wy ] = soft_shrinkage( phiFD.y, abs_phiFD.y, threshr );
                [ phiFD.x, Wx ] = soft_shrinkage( phiFD.x, abs_phiFD.x, threshc );

            end

        else

            error( 'SPECIFY CORRECT TYPE OF THRESHOLDING (''isotropic'' or ''anisotropic'' for sFDxy.nzmask), EXITING )' )

        end

        %=======================

        [ Pi_s_phi( :, :, pp, ss ) ] = iedgedetect_FDxy( phiFD, sFDxy, sol.sz );
        
        

    end
    
end

R_s_phi = 2 * Pi_s_phi - phi;

%==================================================================================================
    
for ss = sol.spos.updateorder
    
    Pi_m_R_s_phi( :, :, :, ss ) = enforce_2DTPAmeas( R_s_phi( :, :, :, ss ), expt.meas.SI( ss ), sol.measLPF, sol );
    
    % Pi_m_phi = enforce_2DTPAmeas( phi, meas, sol.measLPF, sol );
    Pi_m_phi( :, :, :, ss ) = enforce_2DTPAmeas( Pi_s_phi( :, :, :, ss ), expt.meas.SI( ss ), sol.measLPF, sol );
    
end

R_m_R_s_phi = 2 * Pi_m_R_s_phi - R_s_phi;

%==================================================================================================

aalpha = 0.0;
bbeta = 0.5;

phi = 0.5 * aalpha * bbeta * ( R_m_R_s_phi + phi ) + ( 1 - bbeta ) * Pi_m_phi + ( 1 - aalpha )  * Pi_s_phi;
% phi = 0.5 * ( 1 - aalpha ) * ( 1 - bbeta ) * ( R_m_R_s_phi + phi ) + bbeta * Pi_m_phi;

%==================================================================================================

