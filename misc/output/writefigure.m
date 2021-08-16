function writefigure( figID )

%export_fig test3.png -m1.0

tic
export_fig test3.png -r90.0
toc

% tic
% saveas( gcf, [ path_name_svg, '2.jpg' ] )
% toc

% print('PeaksSurface','-dpng')


%print( [ path, name_svg, '.jpg' ], '-djpeg', '-r150' )

% print( [ path, name_svg, '.png' ], '-dpng' )
% print( [ path, name_svg, '.png' ], '-dpng', '-r200' )
% 
% print( [ path, name_svg, '.svg' ], '-dsvg' )

% plot2svg( [ path_name_svg, '.svg' ] ); 
% delete( [ path, name_svg, '001', '.png' ] ); 


