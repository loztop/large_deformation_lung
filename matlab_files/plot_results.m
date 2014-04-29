clear all;
close all;
load matfiles/unconfined_non_p001_726_250.mat

hFig=figure;
set(hFig, 'Position', [100 100 700 400])
set([gca]             , ...
    'FontSize'   , 12           );



num_hand=plot(num_x,num_y,'ko','LineWidth',2,'MarkerSize',7)
hold on;
an_hand=plot(b_x,b_y,'-k','LineWidth',2,'MarkerSize',8)
hold on;
 axis([0 1.2 0.1 0.5]);


%Add the legend and labels
title('Unconfined compression relaxation test','interpreter','latex','FontSize',19);
xlabel('Nondimensional time $(t/t_{g})$ ','interpreter','latex','FontSize',19)
ylabel('Radial displacement $(u/a\epsilon_{0})$ ','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [num_hand,an_hand], ...
  'Numerical solution' ,'Analytical solution',...
  'location', 'NorthEast','interpreter','latex' );

h=hFig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

output_plot_filename='pdf_plots/unconfined_large'
print(h,output_plot_filename,'-dpdf','-r0')

