%clear all;
%close all;


tec='.tec';

prefix='1447TV1.1ke';
Nodes=1447; 
FRC=1.241; %1447
REF=0.9; %1447
T=4;


prefix='5036_0.01W_80NT_8T_';
Nodes=5036; 
FRC=1.477; %1447
REF=0.907; %1447
T=8;


prefix='5036_Cr_0.5_0_';
Nodes=1447; 
FRC=1.477; %1447
REF=0.907; %1447
T=10.4;

NT=[1:1:52];

dt=T/NT(end);

%x_axis_str='Percentage of original Young''s modulus';
x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point/plot_data/';
%out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

%home
%base='/home/loztop/mount_point/mount_point/plot_data/';
%out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';






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

%%Stress
sige(i,:)=data(:,32)';


i
end

mean_s_p=mean(s_p');
std_s_p=std(s_p');

mean_J=mean(J');
std_J=std(J');

mean_sig=mean(sig');
std_sig=std(sig');

mean_sige=mean(sige');
std_sige=std(sige');

t=NT.*dt;

t=[t];


%%calculate volume
Vol=(mean_J-mean_J(1)).*REF;
Vol=(mean_J-1).*REF;

%%calculate flow rate
Frate=(Vol(1:end)-[0,Vol(1:end-1)])./dt


%plotting vars
axisF=23;
fontF=22;

%close all;

hFig=figure;
subplot(2,3,1)
plot(t,[Vol(1:end)],'k-','LineWidth',2,'MarkerSize',8);
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Volume (L)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
% set(gca,'XTick', 0:0.5:4);
%set(gcf, 'PaperSize', [16 15])
%set(gcf, 'PaperPosition', [-1.3 0.2 18 15]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='volume_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')




%hFig=figure;
subplot(2,3,2)
plot(t,Frate,'k-','LineWidth',2,'MarkerSize',8);
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Flow rate (L/s)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
 %set(gca,'XTick', 0:0.5:4);
%set(gcf, 'PaperSize', [16 15])
%set(gcf, 'PaperPosition', [-2 0.2 19  15.5]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='flowrate_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')



 

%hFig=figure;
subplot(2,3,3)
plot(mean_sig(40:end),Vol(40:end),'k-','LineWidth',2,'MarkerSize',8);
hold all;
%arrowhead([mean_sig(9) mean_sig(10)],[Vol(9) Vol(10)],'k',[],2)
hold all;
%arrowhead([mean_sig(29) mean_sig(30)],[Vol(29) Vol(30)],'k',[],2)

set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Volume (L)','interpreter','latex','FontSize',axisF);
xlabel('Mean elastic recoil (Pa)','interpreter','latex','FontSize',axisF);
% set(gca,'XTick', 0:500:3000);
%  set(gca,'YTick', 0:0.1:0.4);
%  ylim([-0.02 0.4]);
%set(gcf, 'PaperSize', [16 15])
%set(gcf, 'PaperPosition', [-1.2 0.2 17.5 15]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='recoil_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')

%close all
%hFig=figure;
subplot(2,3,4)

plot(t,[mean_s_p],'k-','LineWidth',2,'MarkerSize',8);
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Mean pressure drop (Pa)','interpreter','latex','FontSize',axisF);
xlabel('time (s)','interpreter','latex','FontSize',axisF);
% set(gca,'XTick', 0:0.5:4);
set(gcf, 'PaperSize', [16 15])
set(gcf, 'PaperPosition', [-3.7 0.2 21.5 15]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='pressure_tracheaz'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')




%hFig=figure;
subplot(2,3,5)

plot(mean_sige,Vol,'k-','LineWidth',2,'MarkerSize',8);
hold all;
%arrowhead([mean_sig(9) mean_sig(10)],[Vol(9) Vol(10)],'k',[],2)
hold all;
%arrowhead([mean_sig(29) mean_sig(30)],[Vol(29) Vol(30)],'k',[],2)

set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Volume (L)','interpreter','latex','FontSize',axisF);
xlabel('Elastic stress (Pa)','interpreter','latex','FontSize',axisF);











figure;
plot(mean_sig(1:end),Vol(1:end),'k-','LineWidth',2,'MarkerSize',8);


figure;
plot(mean_s_p(1:end),Vol(1:end),'k-','LineWidth',2,'MarkerSize',8);










