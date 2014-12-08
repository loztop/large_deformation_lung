clear all;
close all;


tec='.tec';


prefix='plot_data/5036_Cr_1_0__4';
Nodes=5036; 
FRC=1.477; %1447
REF=0.907; %1447
 
prefix_res='data/5036_H_air_res_0';


NT=[0:1:20];


%work
base='/auto/users/lorenzb/mount_point/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

%home
%base='/home/loztop/mount_point/mount_point/';
%out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';



 i=1;
 
fnameC=strcat(base,prefix,tec)
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


%get airway_tree resistance
fnameC=strcat(base,prefix_res,tec)
importfile(fnameC);
header_info=textdata{3};
TEMP=data;
clear data;
data_res=TEMP(1:Nodes,:);
R=data_res(:,31)';


i

JF=1.3718;
J=((J-JF)./JF)+1;


mean_s_p=mean(s_p');
std_s_p=std(s_p');

mean_J=mean(J');
std_J=std(J');

mean_sig=mean(sig');
std_sig=std(sig');

mean_sige=mean(sige');
std_sige=std(sige');



%%calculate volume
Vol=(mean_J-mean_J(1)).*REF;

close all;


%remove J=1.163
Rscat=R;
Pscat=s_p;

Jscat=J;
idx=find((J<1.1634)&(J>1.1632));
Jscat(idx)=[];
Rscat(idx)=[];
Pscat(idx)=[];


%%calculate correlation coefficients

[r,p]=corrcoef(Rscat',Jscat')
[r,p]=corrcoef(Rscat',Pscat')

corr(Rscat',Jscat')

x=[1:500];
X=x'+600*rand(500,1);
Y=0.5*x'+400*rand(500,1)
[r,p]=corrcoef(X,Y)
close all;
figure; 
scatter(X,Y);

%Plotting business
%plotting vars
axisF=23;
fontF=22;

hFig=figure;
scatter(Rscat,Jscat,'k');
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('${J}_{V}$','interpreter','latex','FontSize',axisF);

xlabel('Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)','interpreter','latex','FontSize',axisF);
 %set(gca,'XTick', 0:0.5:4);
  ylim([1.125 1.184]);
 set(gca,'YTick', 1.13:0.01:1.18);
set(gcf, 'PaperSize', [19 15.4])
set(gcf, 'PaperPosition', [0.2 0.5 19 16]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='vent_res_scatter'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')



hFig=figure;
scatter(Rscat,Pscat,'k');
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Pressure ($\mbox{Pa}$)','interpreter','latex','FontSize',axisF);
xlabel('Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)','interpreter','latex','FontSize',axisF);
 %set(gca,'XTick', 0:0.5:4);
  ylim([-230 20]);
 set(gca,'YTick', -200:50:0);
set(gcf, 'PaperSize', [19 15.4])
set(gcf, 'PaperPosition', [-0.5 0.5 19.5 16]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='pressure_res_scatter'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')


hFig=figure;
scatter(Jscat,Pscat,'k');
set([gca]             , ...
    'FontSize'   , fontF           );
ylabel('Pressure ($\mbox{Pa}$)','interpreter','latex','FontSize',axisF);
xlabel('Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)','interpreter','latex','FontSize',axisF);
 %set(gca,'XTick', 0:0.5:4);
  ylim([-230 20]);
 set(gca,'YTick', -200:50:0);
set(gcf, 'PaperSize', [19 15.4])
set(gcf, 'PaperPosition', [-0.5 0.5 19.5 16]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='pressure_res_scatter'
%print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')



