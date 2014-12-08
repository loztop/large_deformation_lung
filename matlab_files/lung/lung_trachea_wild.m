%clear all;
%close all;

tec='.tec';

prefix='5036_H_80NT_8T_';

prefix='1447_H_NT40_8T_2cyc';
%x_axis_str='Percentage of original Young''s modulus';
x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point/plot_data/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

%home
base='/home/loztop/mount_point/mount_point/plot_data/';
out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

prefix='1447_test';
prefix='5036_pv_8T_2cyc_80NT';

%Nodes=1447;  %1447
%FRC=1.24;  %1447
%REF=9.049;  %1447

Nodes=5036;  
FRC=1.477;  
REF=1.07;  




T=4;


NT=[0:1:80];

for i=1:length(NT)
i
%Load time step data    
fnameC=strcat(base,prefix,'_',num2str(NT(i)),tec)

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


%Pressure
s_p(i,:)=data(:,7)';

%%Jacobian
J(i,:)=data(:,31)';

%%Stress
sig(i,:)=data(:,33)';

i
end



mean_s_p=mean(s_p');
std_s_p=std(s_p');

mean_J=mean(J');
std_J=std(J');

mean_sig=mean(sig');
std_sig=std(sig');

dt=8/NT(end);


TF=80;
PS=41;
PE=81;

t=NT.*dt;

t=[t];


%%calculate volume
%Vol=(mean_J-mean_J(1)).*REF;
Vol=2*[(mean_J).*REF - FRC];
Vol=2*[(mean_J-mean_J(1)).*REF];

%%calculate flow rate
Frate=(Vol(1:end)-[0,Vol(1:end-1)])./dt

%plotting vars
axisF=23;
fontF=22;

close all;

hFig=figure;
plot(t(PS:PE),Vol(PS:PE),'k-','LineWidth',2,'MarkerSize',8);

%plot(t(1:21),Vol(1:21),'k-','LineWidth',2,'MarkerSize',8);
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Volume (L)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
 set(gca,'XTick', 4:1:8);
  set(gca,'YTick', 0:0.1:0.6);
   ylim([-0.02 0.6]);
set(gcf, 'PaperSize', [16 15])
set(gcf, 'PaperPosition', [-1.3 0.2 18 15]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='volume_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')

hFig=figure;
plot(t(PS:PE),Frate(PS:PE),'k-','LineWidth',2,'MarkerSize',8);
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Flow rate (L/s)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
 set(gca,'XTick', 4:1:8);
  set(gca,'YTick', -0.5:0.25:0.5);
set(gcf, 'PaperSize', [16 15])
set(gcf, 'PaperPosition', [-2 0.2 19  15.5]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='flowrate_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')


hFig=figure;
plot(t(PS:PE),[mean_s_p(PS:PE)],'k-','LineWidth',2,'MarkerSize',8);

set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Mean pressure drop (Pa)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
 set(gca,'XTick', 4:1:8);
 
 set(gca,'YTick', -50:25:50);

set(gcf, 'PaperSize', [16 15])
set(gcf, 'PaperPosition', [-3.7 0.2 21.5 15]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='pressure_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')


