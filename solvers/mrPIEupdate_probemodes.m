function [ P ] = mrPIEupdate_probemodes( P, T, phi, sol )

% abs2_T = abs( T ) .^ 2;
% max_abs_sample_abs2 = max( abs2_T( : ));

max_abs_sample_abs2 = max( abs( T( : )) .^ 2 );

if ~isfield( sol.spos, 'updateorder' )
    
    updateorder = randperm( length( sol.spos.indxsubset ), round( 1.0 * length( sol.spos.indxsubset )));
else
    
    updateorder = sol.spos.updateorder;
    
end

%==================================================================================================

tt = 0;
Tmom = 10;
eta_a = 0.7;
eta_b = 0.7;
gamma_m = 0.3;
P_j_plus_1_minusT = P;
nu_j_minus_T = zeros( [ sol.sz.sz, sol.probe.scpm.N ], 'single');

for ss = updateorder       % order * DOES * matter here !!!

%     rs = +1 * sol.spos.rs( ss, : );                 % !!!!!!!!!!!

    %========================================
        
% PRINT OUT % DONE, EVERY 10%

    %========================================
       
    TFview = getview_2DsampleTF( T, sol.sample.vs, sol.spos.rs( ss, : ), sol.spos.shifttype );

%     P = P + conj( TFview ) .* ( phi( :, :, :, ss ) - P .* TFview ) / max_abs_sample_abs2;
    
    aa = 0.9;
%     P = P + gamma_m * conj( TFview ) .* ( phi( :, :, :, ss ) - P .* TFview ) ./ ( aa * max_abs_sample_abs2 + ( 1 - aa) * abs( TFview ) .^ 2 );

    P = P + gamma_m * conj( TFview ) .* ( phi( :, :, :, ss ) - P .* TFview ) ./ max_abs_sample_abs2;
    
    %========================================
    
    P_prime = P;
    tt = tt + 1;

    if ( mod( tt, Tmom ) == 0 )
        
        
       
        nu_j = eta_a * nu_j_minus_T + P_prime - P_j_plus_1_minusT;

%         P = P + nu_j;
        P = P_prime + eta_b * nu_j;
         

        
        tt = 0;
        nu_j_minus_T = nu_j;
        P_j_plus_1_minusT = P;
        
        

    end
    
%     for pp = 1 : sol.probe.scpm.N
% 
%         update_term_p = conj( TFview ) .* ( phi( :, :, pp, ss ) - P( :, :, pp ) .* TFview );
% 
%         P( :, :, pp ) = P( :, :, pp ) + 1.0 * update_term_p / max_abs_sample_abs2;
% 
%     end
    

    
end

5;