clear all;


hFig=figure;
set(hFig, 'Position', [100 100 900 700])
set([gca]             , ...
    'FontSize'   , 12           );



load unconfined_0_001.mat
num_hand_0_001d=plot(num_x,num_y,'ko','LineWidth',2,'MarkerSize',6)
hold on;
an_hand=plot(b_x,b_y,'-k','LineWidth',2,'MarkerSize',8)
hold on;
 axis([0 1.2 0.1 0.5]);
hold all


load unconfined_0_1.mat
num_hand_0_1d=plot(num_x,num_y,'r--','LineWidth',2,'MarkerSize',7)
hold on;
an_hand=plot(b_x,b_y,'-k','LineWidth',2,'MarkerSize',8)
hold on;
hold all

load unconfined_1.mat
num_hand_1d=plot(num_x,num_y,'b.-','LineWidth',2,'MarkerSize',7)
hold on;
an_hand=plot(b_x,b_y,'-k','LineWidth',2,'MarkerSize',8)
hold on;
hold all

%Add the legend and labels
title('Unconfined compression relaxation test','interpreter','latex','FontSize',19);
xlabel('Nondimensional time $(t/t_{g})$ ','interpreter','latex','FontSize',19)
ylabel('Radial displacement $(u/a\epsilon_{0})$ ','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [an_hand,num_hand_0_001d,num_hand_0_1d,num_hand_1d], ...
  'Analytical solution  ','Numerical solution, \delta=0.001  .' ,'Numerical solution, \delta=0.1  ' ,'Numerical solution, \delta=1  ' ,...
  'FontSize',19,'location', 'NorthEast' );

h=hFig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%output_plot_filename='~/Dropbox/Dphil/linear_poro_paper/diagrams/unconfined_results_delta'
%print(h,output_plot_filename,'-dpdf','-r0')

