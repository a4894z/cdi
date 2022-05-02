function hsv_img = complex2hsv( plot_me, log_absrptn, globl_phs, unwrap )

% Takes a complex valued array, and represents it in an HSV colorspace,
% i.e. phase --> hue and magnitude --> value. The saturation component
% is set to 1. 
%
% plot_me   -->   a complex valued 2D array.
%
% globl_phs -->   a real valued scalar, used to give plot_me a global 
%                 phase shift, which in turn will change the hue of the 
%                 image.
%
% log_absrptn --> when we defined the transmission function using some 
%                 image, we either used the value component of the image
%                 or exp(-1 * value component of the image). 
%                 The exp(-1 ... ) is also what we would realistically
%                 be recovering in an experiment. So, if you set this
%                 value to a real valued scalar of 1, you'll take the 
%                 log of plot_me before setting it to the value component.
%
%
% Outputs:
% hsv_img   -->   the complex valued array in an HSV colorspace. Dimenensions
%                 are Ny x Nx x 3.
%
%


% hsv_img( :, :, 1 ) = ( angle( plot_me ) + pi ) / ( 2 * pi );
% hsv_img( :, :, 2 ) = ones( size( plot_me ), 'single' );

% abs_plot_me = abs( plot_me );
% hsv_img( :, :, 3 ) = modulus_limits_project( abs( plot_me ), [ 0.0, 1.0 ] );
% hsv_img( :, :, 3 ) = abs_plot_me / ( 1e-7 + max( abs_plot_me( : )));


%========

if ( unwrap == true )
    
    temp1 = phase_unwrap( angle( plot_me * exp( 1i * 2 * pi * globl_phs )));
    temp1 = temp1 - min( temp1( : ));
    hsv_img( :, :, 1 ) = temp1 / max( temp1( : ));
    
else

    hsv_img( :, :, 1 ) = ( angle( plot_me * exp( 1i * 2 * pi * globl_phs )) + pi ) / ( 2 * pi );
    
end

%========

Sz = size( plot_me );
hsv_img( :, :, 2 ) = ones( Sz( 1 ), Sz( 2 ));

%========
abs_plot_me = abs( plot_me );

if log_absrptn == 1
  
  temp4 = -1 * log( abs( plot_me )); 
  temp4( isinf( temp4 )) = 0; 
  hsv_img( :, :, 3 ) = temp4 / ( 1e-7 + max( temp4( : )));
%   hsv_img( :, :, 3 ) = temp4;
  
else
  
%   hsv_img( :, :, 3 ) = modulus_limits_project( abs_plot_me, [ 0.0, 1.0 ] );
  hsv_img( :, :, 3 ) = abs_plot_me / ( 1e-7 + max( abs_plot_me( : )));
%   hsv_img( :, :, 3 ) = abs_plot_me;
  
end

