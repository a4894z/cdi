function [ Z_bounded ] = array_value_limits( Z, lo_lim, hi_lim )

% RENAME TO ARRAY LIMITS

%Z = angle( x ) + pi; % now phase is between 0 and 2 * pi instead of -pi and +pi


Z_less_than_hi = ( Z < hi_lim );
Z_grtr_than_lo = ( Z > lo_lim );

Z_bounded =  ( Z_less_than_hi & Z_grtr_than_lo ) .* Z + not( Z_less_than_hi ) * hi_lim + not( Z_grtr_than_lo ) * lo_lim;
%Z_bounded =  (Z_less_than_hi & Z_grtr_than_lo) .* Z;   
    
    

% % case 2, swap lo_lim to hi_lim, vice versa
% if (sign( lo_lim ) > 0) && (sign( hi_lim ) <= 0)
%     temp1 = lo_lim; lo_lim = hi_lim; hi_lim = temp1; 
% end
% 
% 
% 
% 
% 
% if (sign( lo_lim ) <= 0) && (sign( hi_lim ) > 0) % case 1
% 
% %     Z_pos = ( Z > 0 );
% %     Z_neg = not( Z_pos );
%     
%     Z_less_than_hi = ( Z <= hi_lim );
%     Z_grtr_than_lo = ( Z >= lo_lim );
%     
% %     Z_bounded = Z .* ( Z_less_than_hi & Z_pos | Z_grtr_than_lo & Z_neg ) + ...
% %         ( not(Z_less_than_hi) .* Z_pos ) * hi_lim + ...
% %         ( not(Z_grtr_than_lo) .* Z_neg ) * lo_lim;
%     
%     Z_bounded = Z .* ( Z_less_than_hi & Z_grtr_than_lo ) + ...
%         not(Z_less_than_hi) * hi_lim + not(Z_grtr_than_lo) * lo_lim;
%     
% %     norm(Z_bounded - Z_bounded2,'fro')
% %     5;
%     
% elseif (sign( lo_lim ) > 0) && (sign( hi_lim ) <= 0) 
%     
% 
% elseif (sign( lo_lim ) <= 0) && (sign( hi_lim ) < 0) % case 3
%         
%         
% elseif (sign( lo_lim ) >= 0) && (sign( hi_lim ) > 0) % case 4
%     
% 
%     Z_less_than_hi = ( Z <= hi_lim );
%     Z_grtr_than_lo = ( Z >= lo_lim );
%     
%     Z_bounded =  (Z_less_than_hi & Z_grtr_than_lo) .* Z + not(Z_less_than_hi) * hi_lim + not(Z_grtr_than_lo) * lo_lim;
%           
% end




%constrain phase:
%{

hi_lim = 1.0*pi/1.0;
lo_lim = -1.0*pi/2000;

temp2 = angle(x);


temp2 = temp2 .* ((temp2<=hi_lim) .* (temp2>=0)) + ...
        temp2 .* ((temp2>=lo_lim) .* (temp2<=0)) + ...
        ((temp2>hi_lim) .* (temp2>=0))*hi_lim + ...
        ((temp2<lo_lim) .* (temp2<=0))*lo_lim;


x_bounded = abs(x) .* exp(1i*temp2);

%}