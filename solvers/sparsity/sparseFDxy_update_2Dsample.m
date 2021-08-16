function [ sample ] = sparseFDxy_update_2Dsample( sample, sparse )
    
    sparse.err_message = '!!!!!!!!!!!!!! sparseFDxy_update_2Dsample !!!!!!!!!!!!!!';

    %==============
    
    switch sparse.threshname

        case 'aniso_abs'

            [ sample ] = sparseFDxy_2Dsample_aniso_abs( sample, sparse );         % sparseFDxy_2Dsample_anisotroptic vs sparseFDxy_2Dsample_anisotroptic_abs

        case 'aniso_phs'

            [ sample ] = sparseFDxy_2Dsample_aniso_phase( sample, sparse );

        case 'iso_abs'

            [ sample ] = sparseFDxy_2Dsample_iso_abs( sample, sparse );

        case 'iso_phs'

            [ sample ] = sparseFDxy_2Dsample_iso_phase( sample, sparse );
            
        case 'aniso_abs_aniso_phs'
         
        [ sample ] = sparseFDxy_2Dsample_aniso_abs( sample, sparse );
        [ sample ] = sparseFDxy_2Dsample_aniso_phase( sample, sparse );
        
%         [ sampleTFabs ] = sparseFDxy_2Dsample_iso_abs( sample, sparse );
%         [ sampleTFphs ] = sparseFDxy_2Dsample_iso_phase( sample, sparse );
%         sample = sampleTFabs .* exp( 1i * sampleTFphs );  
            
        case 'iso_abs_iso_phs'

%             [ sampleTFabs ] = sparseFDxy_2Dsample_iso_abs( sample, sparse );
%             [ sampleTFphs ] = sparseFDxy_2Dsample_iso_phase( sample, sparse );
%             sample = sampleTFabs .* exp( 1i * sampleTFphs );

            [ sample ] = sparseFDxy_2Dsample_iso_abs( sample, sparse );
            [ sample ] = sparseFDxy_2Dsample_iso_phase( sample, sparse );
            
        case 'aniso_re_aniso_im'

            [ sample ] = sparseFDxy_2Dsample_aniso_re_aniso_im( sample, sparse );

        case 'iso_re_iso_im'

            [ sample ] = sparseFDxy_2Dsample_iso_re_iso_im( sample, sparse );

        otherwise
 
            warning('Invalid thresholding type, none performed.')

    end

end

%==================================================================================================
  
function [ sample ] = sparseFDxy_2Dsample_aniso_abs( sample, sparse )

    [ edges_sampleTF ] = edgedetect_FDxy( sample, sparse.s2DFDxy );

    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    abs_edges_sampleTF.y = abs( edges_sampleTF.y );
    abs_edges_sampleTF.x = abs( edges_sampleTF.x );

    lamr = find_thresh_from_sparsitylevel( abs_edges_sampleTF.y, sparse.lvl );
    lamc = find_thresh_from_sparsitylevel( abs_edges_sampleTF.x, sparse.lvl );

    switch sparse.threshtype
        
        case 'h'
    
            [ edges_sampleTF.y, Wr ] = hard_shrinkage( edges_sampleTF.y, abs_edges_sampleTF.y, lamr );
            [ edges_sampleTF.x, Wc ] = hard_shrinkage( edges_sampleTF.x, abs_edges_sampleTF.x, lamc );
    
        case 's'
            
            [ edges_sampleTF.y, Wr ] = soft_shrinkage( edges_sampleTF.y, abs_edges_sampleTF.y, lamr );
            [ edges_sampleTF.x, Wc ] = soft_shrinkage( edges_sampleTF.x, abs_edges_sampleTF.x, lamc );

        otherwise

            error( sparse.err_message )

    end
    
    [ sample ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );

end

%==================================================================================================

function [ sample ] = sparseFDxy_2Dsample_aniso_phase( sample, sparse )

    [ edges_sampleTF ] = edgedetect_FDxy( angle( sample ), sparse.s2DFDxy );

    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    abs_edges_sampleTF.y = abs( edges_sampleTF.y );
    abs_edges_sampleTF.x = abs( edges_sampleTF.x );

    lamr = find_thresh_from_sparsitylevel( abs_edges_sampleTF.y, sparse.lvl );
    lamc = find_thresh_from_sparsitylevel( abs_edges_sampleTF.x, sparse.lvl );

    switch sparse.threshtype
        
        case 'h'
    
            [ edges_sampleTF.y, Wr ] = hard_shrinkage( edges_sampleTF.y, abs_edges_sampleTF.y, lamr );
            [ edges_sampleTF.x, Wc ] = hard_shrinkage( edges_sampleTF.x, abs_edges_sampleTF.x, lamc );
    
        case 's'
            
            [ edges_sampleTF.y, Wr ] = soft_shrinkage( edges_sampleTF.y, abs_edges_sampleTF.y, lamr );
            [ edges_sampleTF.x, Wc ] = soft_shrinkage( edges_sampleTF.x, abs_edges_sampleTF.x, lamc );

        otherwise

            error( sparse.err_message )

    end
    
    [ tmp0 ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );

    sample = abs( sample ) .* exp( 1i * tmp0 );
      
end
        
%==================================================================================================
        
function [ sample ] = sparseFDxy_2Dsample_iso_abs( sample, sparse )

    [ edges_sampleTF ] = edgedetect_FDxy( sample, sparse.s2DFDxy );

    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    abs_edges_sampleTF.y = abs( edges_sampleTF.y );
    abs_edges_sampleTF.x = abs( edges_sampleTF.x );

    iso_abs_edges = sqrt( abs_edges_sampleTF.y .^ 2 + abs_edges_sampleTF.x .^ 2 );

    lam = find_thresh_from_sparsitylevel( iso_abs_edges, sparse.lvl );
    
    switch sparse.threshtype
        
        case 'h'
    
            [ edges_sampleTF.y, Wr ] = hard_shrinkage( edges_sampleTF.y, iso_abs_edges, lam );
            [ edges_sampleTF.x, Wc ] = hard_shrinkage( edges_sampleTF.x, iso_abs_edges, lam );

        case 's'
            
            [ edges_sampleTF.y, Wr ] = soft_shrinkage( edges_sampleTF.y, iso_abs_edges, lam );
            [ edges_sampleTF.x, Wc ] = soft_shrinkage( edges_sampleTF.x, iso_abs_edges, lam );

        otherwise

            error( sparse.err_message )

    end

    [ sample ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );

end

%==================================================================================================

function [ sample ] = sparseFDxy_2Dsample_iso_phase( sample, sparse )


    [ edges_sampleTF ] = edgedetect_FDxy( angle( sample ), sparse.s2DFDxy );

    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    abs_edges_sampleTF.y = abs( edges_sampleTF.y );
    abs_edges_sampleTF.x = abs( edges_sampleTF.x );

    iso_angle_edges = sqrt( abs_edges_sampleTF.y .^ 2 + abs_edges_sampleTF.x .^ 2 );

    lam = find_thresh_from_sparsitylevel( iso_angle_edges, sparse.lvl );
    
    switch sparse.threshtype
        
        case 'h'
    
            [ edges_sampleTF.y, Wr ] = hard_shrinkage( edges_sampleTF.y, iso_angle_edges, lam );
            [ edges_sampleTF.x, Wc ] = hard_shrinkage( edges_sampleTF.x, iso_angle_edges, lam );

        case 's'
            
            [ edges_sampleTF.y, Wr ] = soft_shrinkage( edges_sampleTF.y, iso_angle_edges, lam );
            [ edges_sampleTF.x, Wc ] = soft_shrinkage( edges_sampleTF.x, iso_angle_edges, lam );

        otherwise

            error( sparse.err_message )

    end
    
    [ tmp0 ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );

    sample = abs( sample ) .* exp( 1i * tmp0 );
    
end
        
%==================================================================================================

function [ sample ] = sparseFDxy_2Dsample_aniso_re_aniso_im( sample, sparse )
  
    [ edges_sampleTF ] = edgedetect_FDxy( sample, sparse.s2DFDxy );

    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    re_edges_sampleTF.r = real( edges_sampleTF.y );
    re_edges_sampleTF.c = real( edges_sampleTF.x );
    im_edges_sampleTF.r = imag( edges_sampleTF.y );
    im_edges_sampleTF.c = imag( edges_sampleTF.x );

    abs_re_edges_sampleTF.r = abs( re_edges_sampleTF.r );
    abs_re_edges_sampleTF.c = abs( re_edges_sampleTF.c );
    abs_im_edges_sampleTF.r = abs( im_edges_sampleTF.r );
    abs_im_edges_sampleTF.c = abs( im_edges_sampleTF.c );

    lamrer = find_thresh_from_sparsitylevel( abs_re_edges_sampleTF.r, sparse.lvl );
    lamrec = find_thresh_from_sparsitylevel( abs_re_edges_sampleTF.c, sparse.lvl );
    lamimr = find_thresh_from_sparsitylevel( abs_im_edges_sampleTF.r, sparse.lvl );
    lamimc = find_thresh_from_sparsitylevel( abs_im_edges_sampleTF.c, sparse.lvl );

    switch sparse.threshtype
        
        case 'h'
    
            [ re_edges_sampleTF.r, Wre.r ] = hard_shrinkage( re_edges_sampleTF.r, abs_re_edges_sampleTF.r, lamrer );
            [ re_edges_sampleTF.c, Wre.c ] = hard_shrinkage( re_edges_sampleTF.c, abs_re_edges_sampleTF.c, lamrec );
            [ im_edges_sampleTF.r, Wim.r ] = hard_shrinkage( im_edges_sampleTF.r, abs_im_edges_sampleTF.r, lamimr );
            [ im_edges_sampleTF.c, Wim.c ] = hard_shrinkage( im_edges_sampleTF.c, abs_im_edges_sampleTF.c, lamimc );

        case 's'
            
            [ re_edges_sampleTF.r, Wre.r ] = soft_shrinkage( re_edges_sampleTF.r, abs_re_edges_sampleTF.r, lamrer );
            [ re_edges_sampleTF.c, Wre.c ] = soft_shrinkage( re_edges_sampleTF.c, abs_re_edges_sampleTF.c, lamrec );
            [ im_edges_sampleTF.r, Wim.r ] = soft_shrinkage( im_edges_sampleTF.r, abs_im_edges_sampleTF.r, lamimr );
            [ im_edges_sampleTF.c, Wim.c ] = soft_shrinkage( im_edges_sampleTF.c, abs_im_edges_sampleTF.c, lamimc );

        otherwise

            error( sparse.err_message )

    end
    
    edges_sampleTF.y = re_edges_sampleTF.r + 1i * im_edges_sampleTF.r;
    edges_sampleTF.x = re_edges_sampleTF.c + 1i * im_edges_sampleTF.c;

    [ sample ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );     

end

%==================================================================================================

function [ sample ] = sparseFDxy_2Dsample_iso_re_iso_im( sample, sparse )
  

    [ edges_sampleTF ] = edgedetect_FDxy( sample, sparse.s2DFDxy );


    edges_sampleTF.y = edges_sampleTF.y .* sparse.support;
    edges_sampleTF.x = edges_sampleTF.x .* sparse.support;

    re_edges_sampleTF.r = real( edges_sampleTF.y );
    re_edges_sampleTF.c = real( edges_sampleTF.x );
    im_edges_sampleTF.r = imag( edges_sampleTF.y );
    im_edges_sampleTF.c = imag( edges_sampleTF.x );

    abs_re_edges_sampleTF.r = abs( re_edges_sampleTF.r );
    abs_re_edges_sampleTF.c = abs( re_edges_sampleTF.c );
    tmp0re = sqrt( abs_re_edges_sampleTF.r .^ 2 + abs_re_edges_sampleTF.c .^ 2 );

    abs_im_edges_sampleTF.r = abs( im_edges_sampleTF.r );
    abs_im_edges_sampleTF.c = abs( im_edges_sampleTF.c );
    tmp0im = sqrt( abs_im_edges_sampleTF.r .^ 2 + abs_im_edges_sampleTF.c .^ 2 );

    lamre = find_thresh_from_sparsitylevel( tmp0re, sparse.lvl );
    lamim = find_thresh_from_sparsitylevel( tmp0im, sparse.lvl );

    switch sparse.threshtype
        
        case 'h'
    
            [ re_edges_sampleTF.r, Wre.r ] = hard_shrinkage( re_edges_sampleTF.r, tmp0re, lamre );
            [ re_edges_sampleTF.c, Wre.c ] = hard_shrinkage( re_edges_sampleTF.c, tmp0re, lamre );
            [ im_edges_sampleTF.r, Wim.r ] = hard_shrinkage( im_edges_sampleTF.r, tmp0im, lamim );
            [ im_edges_sampleTF.c, Wim.c ] = hard_shrinkage( im_edges_sampleTF.c, tmp0im, lamim );

        case 's'
            
            [ re_edges_sampleTF.r, Wre.r ] = soft_shrinkage( re_edges_sampleTF.r, tmp0re, lamre );
            [ re_edges_sampleTF.c, Wre.c ] = soft_shrinkage( re_edges_sampleTF.c, tmp0re, lamre );
            [ im_edges_sampleTF.r, Wim.r ] = soft_shrinkage( im_edges_sampleTF.r, tmp0im, lamim );
            [ im_edges_sampleTF.c, Wim.c ] = soft_shrinkage( im_edges_sampleTF.c, tmp0im, lamim );

        otherwise

            error( sparse.err_message )

    end

    edges_sampleTF.y = re_edges_sampleTF.r + 1i * im_edges_sampleTF.r;
    edges_sampleTF.x = re_edges_sampleTF.c + 1i * im_edges_sampleTF.c;

    [ sample ] = iedgedetect_FDxy( edges_sampleTF, sparse.s2DFDxy );            

end

