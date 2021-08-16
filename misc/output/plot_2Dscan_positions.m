function plot_2Dscan_positions( sposA, labelindxA, sposB, labelindxB )

%==================================================================================================

hold on 

if ~isempty( sposA ) 
    
    scatter( sposA( :, 2 ), sposA( :, 1 ), 30, 'o', ...
                                'MarkerEdgeColor', [ 0.0, 0.0, 0.0 ], ...
                                'MarkerEdgeAlpha', 0.3, ...
                                'MarkerFaceAlpha', 0.3, ...
                                'LineWidth', 1.5 );

    % scatter( sposA( :, 2 ), sposA( :, 1 ), 30, 'o', ...
    %                             'MarkerEdgeColor', [ 0.0, 0.0, 0.0 ], ...
    %                             'MarkerFaceColor', [ 0.0, 0.0, 0.0 ], ...
    %                             'LineWidth', 1.5 );

end


if ~isempty( labelindxA ) 

    sposA = double( sposA );
    labelindxA = double( labelindxA );
    
    for ii = 1 : 1 : length( labelindxA )



% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% 
%         text( sposA( labelindxA( ii ), 2 ), sposA( labelindxA( ii ), 1 ) + 0e-0 , num2str( labelindxA( ii ), '%d' ), ...
%                                                     'FontSize', 8, ...
%                                                     'FontWeight', 'Bold', ...
%                                                     'Color', [ 0.8, 0.0, 0.0 ], ...
%                                                     'HorizontalAlignment', 'center' );
                                                

        text( sposA( ii, 2 ), sposA( ii, 1 ) + 0.5e-6, num2str( labelindxA( ii ), '%d' ), ...
                                                    'FontSize', 8, ...
                                                    'FontWeight', 'Bold', ...
                                                    'Color', [ 0.8, 0.0, 0.0 ], ...
                                                    'HorizontalAlignment', 'center' );
                                                                                       

    end


end

%==================================================================================================

if ~isempty( sposB )

    scatter( sposB( :, 2 ), sposB( :, 1 ), 10, 'o', ...
                                'MarkerEdgeColor', [ 0.0, 0.7, 0.0 ], ...
                                'MarkerFaceColor', [ 0.0, 0.7, 0.0 ], ...
                                'MarkerEdgeAlpha', 0.6, ...
                                'MarkerFaceAlpha', 0.6, ...
                                'LineWidth', 1.5 );
                            
end




if ~isempty( labelindxB )


    sposB = double( sposB );
    labelindxB = double( labelindxB );

    for ii = 1 : 1 : length( labelindxB )


% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% % OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
% 
%         text( sposB( labelindxB( ii ), 2 ), sposB( labelindxB( ii ), 1 ) + 0.35e-6 , num2str( labelindxB( ii ), '%d' ), ...
%                                                     'FontSize', 8, ...
%                                                     'FontWeight', 'Bold', ...
%                                                     'Color', [ 0.0, 0.0, 0.8 ], ...
%                                                     'HorizontalAlignment', 'center' );   
                                                
                                                
        text( sposB( ii, 2 ), sposB( ii, 1 ) + 0 * 0.2e-6 , num2str( labelindxB( ii ), '%d' ), ...
                                                    'FontSize', 8, ...
                                                    'FontWeight', 'Bold', ...
                                                    'Color', [ 0.0, 0.0, 0.8 ], ...
                                                    'HorizontalAlignment', 'center' );  

                            
                                                
    hold on 
                         
    xl( 1 ) = sposA( labelindxB( ii ), 2 ); yl( 1 ) = sposA( labelindxB( ii ), 1 );
    xl( 2 ) = sposB( ii, 2 ); yl( 2 ) = sposB( ii, 1 );
    line( xl, yl, 'LineWidth', 1, 'Color', [1.0 0.0 0.0] );
    
    hold off                                      

    end

end




% hold off
% 
% if all( size( sposA ) == size( sposB ))
%     
%     hold on 
%                          
%     xl( 1, : ) = sposA( :, 2 ); yl( 1, : ) = sposA( :, 1 );
%     xl( 2, : ) = sposB( :, 2 ); yl( 2, : ) = sposB( :, 1 );
%     line( xl, yl, 'LineWidth', 1, 'Color', [1.0 0.0 0.0] );
%     
%     hold off
%     
% end

% daspect([1 1 1])   
   
%==================================================================================================



%{

x = 1:500; 
y = sind(x); 
plot(x,y,'linewidth',3) 
axis tight; 
hold on
I = imread('peppers.png'); 
h = image(xlim,-ylim,I); 
uistack(h,'bottom')

%}



%{

% This creates the 'background' axes
ha = axes('units','normalized', 'position',[0 0 1 1]);
        
% Move the background axes to the bottom
uistack(ha,'bottom');

% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
I=imread('eight.tif');

hi = imagesc(I);

colormap gray

% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
set(ha,'handlevisibility','off', 'visible','off')
        
% Now we can use the figure, as required.
% For example, we can put a plot in an axes
axes('position',[0.3,0.35,0.4,0.4])

imagesc( imresize( I, 0.5 ))

%}