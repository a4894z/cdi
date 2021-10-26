function demo_shadedErrorBar

%{

cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/shadedErrorBar

%}

%====================================================================================================================================================

% figure
% %Plot the mean and standard deviation then overlay the raw data
% y = randn( 30, 80 ) * 5;
% x = ( 1 : size( y, 2 )) - 40;
% 
% yP = sin( linspace( -2 * pi, 2 * pi, length( x )) ) * 20;
% 
% % y = bsxfun( @plus, y, yP ) + 60;
% y = y + yP + 60;
% 
% shadedErrorBar( x, y, { @mean, @std }); 
% 
% hold on
% plot( x, y, '.', 'color', [ 0.5, 0.5, 0.95 ] )
% 
% hold off
% 
% grid on
% 
% 5;
% 
% %====================================================================================================================================================
% 
% figure
% %Overlay different lines (transparent) lines and change their properties
% hold on
% 
% plot2styles = { '-b'; '-g'; '-r' };
% 
% for i = 1 : 3
%     
%   plt2invis = plot( 0, 0, plot2styles{ i });
%   
%   set( plt2invis, 'visible', 'off' ); 
%   
% end
%   
% x = ( 1 : size( y, 2 )) - 40;
% y = ones( 30, 1 ) * x;
% 
% y = y + 0.06 * y .^ 2 + randn( size( y )) * 10;
% 
% 
% 
% 
% 
% % min_max_bounds = @( n ) bounds( n, [], 1 );
% 
% % figure; 
% 
% shadedErrorBar( x, y, { @mean, @std }, 'lineprops', '-b', 'transparent', true, 'patchSaturation', 0.5 )
% % shadedErrorBar( x, y, { @mean, min_max_bounds }, 'lineprops', '-b', 'transparent', true, 'patchSaturation', 0.5 )
% 
% % figure; 
% % plot(  x, max(y,[],1) ); hold on; plot(  x, min(y,[],1) ); plot(  x, mean(y,1) ); 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% shadedErrorBar( x, 2 * y + 20, { @mean, @std }, 'lineprops', { '-o', 'Color', [0, 0.7, 0], 'MarkerFaceColor', [0, 0.7, 0] }, 'transparent', true, 'patchSaturation', 0.4 );
% 
% y = randn( 30, 80 ) * 5; 
% x = ( 1 : size( y, 2 )) - 40;
% 
% yP = sin( linspace( -2 * pi, 2 * pi, length( x ))) * 20;
% y  = bsxfun( @plus, y, yP ) + 60;
% 
% %Make red line non-transparent
% shadedErrorBar( x, y, { @mean, @std }, 'lineprops', '-r', 'transparent', true, 'patchSaturation', 0.075 )
% hold off
% 
% legend( [ '#1'; '#2'; '#3' ], 'location', 'northwest' )
% legend
% grid on
% 
% 5;

%====================================================================================================================================================

figure
% Post-hoc modifications of line properties


y=randn(30,80)*5; 
x=(1:size(y,2));
yP = sin( linspace(-2*pi,2*pi,length(x)) )*20;
y = bsxfun(@plus,y,yP);


% %Set face and edge properties
% if (sum( size(ver('MATLAB'))) > 0  )
%  s = shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r')
%  set(s.edge,'LineWidth',2,'LineStyle',':')
%  s.mainLine.LineWidth = 5;
%  s.patch.FaceColor = [0.5,0.25,0.25];
% elseif (sum(size(ver('Octave'))) > 0)
%  s = shadedErrorBar(x, y, {@mean,@std}, 'lineprops', {'-r','LineWidth',5,'LineStyle',':'});
% end

s = shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r');
set(s.edge,'LineWidth',2,'LineStyle',':')
s.mainLine.LineWidth = 5;
s.patch.FaceColor = [0.5,0.25,0.25];



hold on

if (sum( size(ver('MATLAB'))) > 0  )
 plot(s.mainLine.XData, s.mainLine.YData,'or','MarkerFaceColor','w')
end

hold off
grid on

set(gca,'XTickLabel',[],'YTickLabel',[])

%====================================================================================================================================================
% Post-hoc modifications of line properties

% figure
% y=randn(256,80)*5; 
% x=(1:size(y,2));
% yP = cos( linspace(-2*pi,2*pi,length(x)) )*10;
% y = bsxfun(@plus,y,yP);
% 
% 
% shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r')
% 
% hold on
% 
% y=mean(y)+16;
% errbar = [2*ones(1,length(x)) ; 4*ones(1,length(x))];
% 
% shadedErrorBar(x, y, errbar, 'lineprops', '-g')
% hold off
