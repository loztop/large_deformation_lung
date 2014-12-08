clear all;
close all;


tec='.tec';

prefix1='1447_pv_2T_2cyc_80NT';
prefix2='1447_pv_8T_2cyc_80NT';
prefix3='1447_pv_32T_2cyc_80NT';
Nodes=1447; 

prefix1='5036_pv_2T_2cyc_80NT';
prefix2='5036_pv_8T_2cyc_80NT';
prefix3='5036_pv_32T_2cyc_80NT';


Nodes=5036; 

%FRC=1.241; %1447
%REF=0.904; %1447

FRC=1.477; %5036
REF=1.077; %5036

T=2;

TF=80;


NT=[0:1:TF];

dt=T/NT(end);

%x_axis_str='Percentage of original Young''s modulus';
x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
%base='/auto/users/lorenzb/mount_point/plot_data/';
%out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

%home
base='/home/loztop/mount_point/mount_point/plot_data/';
out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';






for i=1:length(NT)
i
%Load time step data    
fnameC1=strcat(base,prefix1,'_',num2str(NT(i)),tec)
fnameC2=strcat(base,prefix2,'_',num2str(NT(i)),tec)
fnameC3=strcat(base,prefix3,'_',num2str(NT(i)),tec)

importfile(fnameC1);
header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:Nodes,:);
s_p(i,:)=data(:,7)';
J(i,:)=data(:,31)';
sig(i,:)=data(:,33)';
sige(i,:)=data(:,32)';




importfile(fnameC2);
header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:Nodes,:);
s_p2(i,:)=data(:,7)';
J2(i,:)=data(:,31)';
sig2(i,:)=data(:,33)';
sige2(i,:)=data(:,32)';



importfile(fnameC3);
header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:Nodes,:);
s_p3(i,:)=data(:,7)';
J3(i,:)=data(:,31)';
sig3(i,:)=data(:,33)';
sige3(i,:)=data(:,32)';

i


t=size(data,1);
end

mean_s_p=mean(s_p');
mean_J=mean(J');
mean_sig=mean(sig');
mean_sige=mean(sige');

mean_s_p2=mean(s_p2');
mean_J2=mean(J2');
mean_sig2=mean(sig2');
mean_sige2=mean(sige2');

mean_s_p3=mean(s_p3');
mean_J3=mean(J3');
mean_sig3=mean(sig3');
mean_sige3=mean(sige3');

t=NT.*dt;

t=[t];


%%calculate volume
%Vol=(mean_J-mean_J(1)).*REF;
Vol=2*[(mean_J).*REF - FRC];

%%calculate flow rate
Frate=2*(Vol(1:end)-[0,Vol(1:end-1)])./dt


%plotting vars
axisF=24;
fontF=22;
PS=41;
PE=TF+1;

Asize=[0.65,0.65];
A1=50;
A2=70;

hFig=figure;
pre1=plot(mean_sig(PS:PE),Vol(PS:PE),'k--','LineWidth',2,'MarkerSize',8);
hold all;
pre2=plot(mean_sig2(PS:PE),Vol(PS:PE),'b-','LineWidth',2,'MarkerSize',8);
hold all;
pre3=plot(mean_sig3(PS:PE),Vol(PS:PE),'r-.','LineWidth',2,'MarkerSize',8);
hold all;


arrowhead([mean_sig(A1) mean_sig(A1+1)],[Vol(A1) Vol(A1+1)],'k',Asize,2)
hold all;
arrowhead([mean_sig(A2) mean_sig(A2+1)],[Vol(A2) Vol(A2+1)],'k',Asize,2)
hold all;

arrowhead([mean_sig2(A1) mean_sig2(A1+1)],[Vol(A1) Vol(A1+1)],'k',Asize,2)
hold all;
arrowhead([mean_sig2(A2) mean_sig2(A2+1)],[Vol(A2) Vol(A2+1)],'k',Asize,2)
hold all;

arrowhead([mean_sig3(A1) mean_sig3(A1+1)],[Vol(A1) Vol(A1+1)],'k',Asize,2)
hold all;
arrowhead([mean_sig3(A2) mean_sig3(A2+1)],[Vol(A2) Vol(A2+1)],'k',Asize,2)
hold all;


% set([gca]             , ...
%    'FontSize'   , fontF           );
% ylabel('Volume (L)','interpreter','latex','FontSize',axisF);
% xlabel('Mean elastic recoil (Pa)','interpreter','latex','FontSize',axisF);
%  set(gca,'XTick', 200:200:1400);
%   set(gca,'YTick', 0.6:0.1:1.2);
%   ylim([0.55 1.25]);
% set(gcf, 'PaperSize', [23.9 18.5])
% set(gcf, 'PaperPosition', [-1.2 0.2 26.4 19]); %Position plot at left hand corner with width 5 and height 5.
% 
% hLegend = legend( ...
%    [pre1,pre2,pre3], ...
%   '1s breathing cycle' ,'4s breathing cycle' ,'16s breathing cycle' ,...
%   'FontSize',18,'location', 'SouthEast' );
% 
% set(hLegend,'interpreter','latex');
% 
% output_plot_filename='recoil_tracheaz'
% 
% print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')














