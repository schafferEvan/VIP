function [axis_limits] = DisplayPoints3DPretty(Model, Scene, TransformedModel)
%%=====================================================================
%% $RCSfile: DisplayPoints3D.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $

% Modified 10/11/2019 ESS
%%=====================================================================
figure
set(gcf,'color','k')
pos=get(gcf,'Position');
heightRescale = 1;
widthRescale = 2;
pos(3)=pos(3)*widthRescale;
pos(4)=pos(4)*heightRescale;
set(gcf,'Position',pos)

subplot(1,2,1);
axis_limits = determine_border(Model, Scene);
set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');

plot3(Model(:,1),Model(:,2),Model(:,3),'r.', 'MarkerSize', 8, 'LineWidth',1.5);
hold on;
plot3(Scene(:,1),Scene(:,2),Scene(:,3),'.','Color',[0,.7,.7], 'MarkerSize', 8, 'LineWidth',1.5);
axis_limits = determine_border(Model, Scene);
xlim(axis_limits(1,:));
ylim(axis_limits(2,:));   
zlim(axis_limits(3,:));   
axis equal;
axis off

subplot(1,2,2);
plot3(TransformedModel(:,1),TransformedModel(:,2),TransformedModel(:,3),'r.', 'MarkerSize', 8, 'LineWidth',1.5);
hold on;
plot3(Scene(:,1),Scene(:,2),Scene(:,3),'.','Color',[0,.7,.7], 'MarkerSize', 8, 'LineWidth',1.5);
axis_limits = determine_border(TransformedModel, Scene);
xlim(axis_limits(1,:));
ylim(axis_limits(2,:));   
zlim(axis_limits(3,:));   
axis equal;
axis off

set(gcf,'Color','k')
drawnow