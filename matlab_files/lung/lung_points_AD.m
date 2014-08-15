clear all;
close all;

tec='.tec';

prefix='1447_C_';

%x_axis_str='Percentage of original Young''s modulus';
x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point/plot_data/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Craf_';

%home
%base='/home/loztop/mount_point/mount_point/plot_data/';
%out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/C_';

Nodes=1447 %for 246 mesh

%Tracking points
%P=[-60,-170,-90;
%    -60,-155,-95;
%    -60,-143,-100;
%    -60,-125,-105];

P=[-60,-170,-90;
    -60,-145,-95;
    -60,-137,-90;
    -60,-125,-100];

NP=size(P,1);

AD=[1,0.5,0.25,0.1];
TD=[0,0,0,0,0];


NT=9;

for i=1:length(AD)
i
%Load time step data    
fnameC=strcat(base,prefix,num2str(AD(i)),'_',num2str(TD(i)),'__',num2str(NT),tec)

importfile(fnameC);

header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:Nodes,:);

%This might change depending on what order vraibles are stored

t=size(data,1);

%reference positions
x=data(:,20);
y=data(:,21);
z=data(:,22);

%loop thorough the nodes to find the closest one to the three points of
%interest

near_dist_P=ones(1,4).*999;
near_P_idx(i,:)=zeros(1,4);
for j=1:Nodes
    
    for k=1:NP
        dist_P(k)=(P(k,1)-x(j)).^2 + (P(k,2)-y(j)).^2 + (P(k,3)-z(j)).^2;
        if(dist_P(k)<near_dist_P(k))
            near_dist_P(k)=dist_P(k);
            near_P_idx(i,k)=j;
        end
    end 
    
end

%Pressure
s_p(i,:)=data(near_P_idx(i,:),7)';
%Jacobian
J(i,:)=data(near_P_idx(i,:),31)';
%'p1' - the toatal stress
sig(i,:)=data(near_P_idx(i,:),32)';

end

AD(1)=1;
AD=AD.*100;
xaxis=AD;

hFig=figure;
set(hFig, 'Position', [100 100 550 430])
set([gca]             , ...
    'FontSize'   , 12           );

p1_hand=plot(xaxis,(s_p(:,1)),'k-s','LineWidth',2,'MarkerSize',8);
hold all;
p2_hand=plot(xaxis,s_p(:,2),'b-o','LineWidth',2,'MarkerSize',8);
hold all;
p3_hand=plot(xaxis,s_p(:,3),'-.');
hold all;
p4_hand=plot(xaxis,s_p(:,4),'g-x','LineWidth',2,'MarkerSize',8);
hold all;
set(gca,'XDir','reverse') 
%Add the legend and labels
title('Alveolar pressure','interpreter','latex','FontSize',19);
xlabel('Percentage of original airway radius','interpreter','latex','FontSize',19)
ylabel('Pressure (Pa)','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [p1_hand,p2_hand,p4_hand], ...
  'Black ball (center disease)','Blue ball (inerface)','Green ball (healthy)',...
  'location', 'SouthWest','interpreter','latex' );
set(hLegend,'interpreter','latex')
set(hLegend,'FontSize',13)

set(gcf, 'PaperPosition', [0 0 18 18]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [18 18]);

axis([0 105 -1600 50]);

output_plot_filename='pressure';
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')




hFig=figure;
set(hFig, 'Position', [100 100 550 430])
set([gca]             , ...
    'FontSize'   , 12           );

p1_hand=plot(xaxis,(J(:,1)),'k-s','LineWidth',2,'MarkerSize',8);
hold all;
p2_hand=plot(xaxis,J(:,2),'b-o','LineWidth',2,'MarkerSize',8);
hold all;
p3_hand=plot(xaxis,J(:,3),'-.');
hold all;
p4_hand=plot(xaxis,J(:,4),'g-x','LineWidth',2,'MarkerSize',8);
hold all;
set(gca,'XDir','reverse') 
%Add the legend and labels
title('Jacobian','interpreter','latex','FontSize',19);
xlabel('Percentage of original airway radius','interpreter','latex','FontSize',19)
ylabel('J','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [p1_hand,p2_hand,p4_hand], ...
  'Black ball (center disease)','Blue ball (inerface)','Green ball (healthy)',...
  'location', 'SouthWest','interpreter','latex' );
set(hLegend,'interpreter','latex')
set(hLegend,'FontSize',13)
set(gcf, 'PaperPosition', [0 0 18 18]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [18 18]);
axis([0 105 1.1 1.25]);

output_plot_filename='jacobian';
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')






hFig=figure;
set(hFig, 'Position', [100 100 550 430])
set([gca]             , ...
    'FontSize'   , 12           );

p1_hand=plot(xaxis,(sig(:,1)),'k-s','LineWidth',2,'MarkerSize',8);
hold all;
p2_hand=plot(xaxis,sig(:,2),'b-o','LineWidth',2,'MarkerSize',8);
hold all;
p3_hand=plot(xaxis,sig(:,3),'-.');
hold all;
p4_hand=plot(xaxis,sig(:,4),'g-x','LineWidth',2,'MarkerSize',8);
hold all;
set(gca,'XDir','reverse') 
%Add the legend and labels
title('Stress','interpreter','latex','FontSize',19);
xlabel(x_axis_str,'interpreter','latex','FontSize',19)
ylabel('Stress (Pa)','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [p1_hand,p2_hand,p4_hand], ...
  'Black ball (center disease)','Blue ball (inerface)','Green ball (healthy)',...
  'location', 'NorthWest','interpreter','latex' );
set(hLegend,'interpreter','latex')
set(hLegend,'FontSize',13)
set(gcf, 'PaperPosition', [0 0 18 18]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [18 18]);

output_plot_filename='stress';
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')
axis([0 105 2000 5000]);
