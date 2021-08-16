function imagescHSV( plot_me, pltopts ) 

% Takes a complex valued array, and represents it in an HSV colorspace,
% i.e. phase --> hue and magnitude --> value. The saturation component
% is set to 1. The resultant HSV image is then plotted using imagesc.
%
% Inputs: 
%
% plot_me   -->   a complex valued 2D array.
%
% global_phase -->   a real valued scalar, used to give plot_me a global 
%                 phase shift, which in turn will change the hue of the 
%                 image.
%
% log_magnitude --> when we defined the transmission function using some 
%                 image, we either used the value component of the image
%                 or exp(-1 * value component of the image). 
%                 The exp(-1 ... ) is also what we would realistically
%                 be recovering in an experiment. So, if you set this
%                 value to a real valued scalar of 1, you'll take the 
%                 log of plot_me before setting it to the value component.
%
%
% Outputs:
%
% none
%
%



%--------------------------------------------------------------------------------------------------

if ~exist( 'pltopts', 'var' )
    
    pltopts = struct;

%     pltopts.xaxis              = 1 : size( plot_me, 2 );
%     pltopts.yaxis              = 1 : size( plot_me, 1 );
%     pltopts.global_phase       = 0;
%     pltopts.unwrap             = false;
%     pltopts.log_magnitude      = 0;
%     pltopts.brightness_scaling = 1;

end

    
if ~isfield( pltopts, 'xaxis' ),                pltopts.xaxis              = 1 : size(plot_me,2); end
if ~isfield( pltopts, 'yaxis' ),                pltopts.yaxis              = 1 : size(plot_me,1); end
if ~isfield( pltopts, 'global_phase' ),         pltopts.global_phase       = 0; end
if ~isfield( pltopts, 'log_magnitude' ),        pltopts.log_magnitude      = 0; end
if ~isfield( pltopts, 'brightness_scaling' ),   pltopts.brightness_scaling = 1; end
if ~isfield( pltopts, 'unwrap' ),               pltopts.unwrap             = false; end


%--------------------------------------------------------------------------------------------------

hsv_img = complex2hsv( plot_me, pltopts.log_magnitude, pltopts.global_phase, pltopts.unwrap );

hsv_img( :, :, 3 ) = pltopts.brightness_scaling * hsv_img( :, :, 3 );

rgb_img = hsv2rgb( hsv_img );

%--------------------------------------------------------------------------------------------------

imagesc( pltopts.xaxis, pltopts.yaxis, rgb_img ); 
%daspect( [1 1 1] ); 

%--------------------------------------------------------------------------------------------------


% print('ycomp.ppm','-dppm','-r600')
 
 
 
end

