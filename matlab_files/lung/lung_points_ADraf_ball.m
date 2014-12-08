clear all;
close all;

tec='.tec';

prefix='5036_Cr_';
%AD=[1,0.6,0.5,0.4,0.35,0.3];
%TD=[0,0,0,0,0,0,0];

AD=[1,0.6,0.5,0.4,0.35];
TD=[0,0,0,0,0,0];

x_axis_str='Percentage of original Young''s modulus';
%x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point2/plot_data/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Wrafeo_';

%home
%base='/home/loztop/mount_point/mount_point/plot_data/';
%out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

Nodes=5036 %for 246 mesh

%Tracking points
%P=[-60,-170,-90;
%    -60,-155,-95;
%    -60,-143,-100;
%    -60,-125,-105];

P=[-60,-170,-90;
    -60,-145,-95;
    -60,-137,-90;
    -60,-125,-100];

C=[-57,-144,-90];

NP=size(P,1);



NT=4;

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
 
        dist_C(j)= sqrt((C(1,1)-x(j)).^2 + (C(1,2)-y(j)).^2 + (C(1,3)-z(j)).^2);
  
    
end

%Pick out nodes for relevant bands of points

b1_c=0;
b2_c=0;
b3_c=0;

for j=1:Nodes
 
        if(dist_C(j)<19)
           b1_c=b1_c+1;
           b1_idx(b1_c)=j;
        end
  
        if(dist_C(j)>19 && dist_C(j)<21)
           b2_c=b2_c+1;
           b2_idx(b2_c)=j;
        end
        
        if(dist_C(j)>24 && dist_C(j)<29 )
           b3_c=b3_c+1;
           b3_idx(b3_c)=j;
        end
    
end

%size(b1_idx)

%Pressure
s_p_b1(i,:)=data(b1_idx,7)';
s_p_b2(i,:)=data(b2_idx,7)';
s_p_b3(i,:)=data(b3_idx,7)';

%%Jacobian
J_b1(i,:)=data(b1_idx,31)';
J_b2(i,:)=data(b2_idx,31)';
J_b3(i,:)=data(b3_idx,31)';

%%Stress
sig_b1(i,:)=data(b1_idx,33)';
sig_b2(i,:)=data(b2_idx,33)';
sig_b3(i,:)=data(b3_idx,33)';

%%Elastic Stress
sige_b1(i,:)=data(b1_idx,32)';
sige_b2(i,:)=data(b2_idx,32)';
sige_b3(i,:)=data(b3_idx,32)';


clear b1_idx
clear b2_idx
clear b3_idx

%%Jacobian
%J(i,:)=data(near_P_idx(i,:),31)';
%%'p1' - the toatal stress
%sig(i,:)=data(near_P_idx(i,:),32)';
i
end

JF=1.3718;
J_b1=((J_b1-JF)./JF)+1;
J_b2=((J_b2-JF)./JF)+1;
J_b3=((J_b3-JF)./JF)+1;


mean_J1=mean(J_b1');
std_J1=std(J_b1');
mean_J2=mean(J_b2');
std_J2=std(J_b2');
mean_J3=mean(J_b3');
std_J3=std(J_b3');

mean_s_p1=mean(s_p_b1');
std_s_p1=std(s_p_b1');
mean_s_p2=mean(s_p_b2');
std_s_p2=std(s_p_b2');
mean_s_p3=mean(s_p_b3');
std_s_p3=std(s_p_b3');

mean_sig1=mean(sig_b1');
std_sig1=std(sig_b1');
mean_sig2=mean(sig_b2');
std_sig2=std(sig_b2');
mean_sig3=mean(sig_b3');
std_sig3=std(sig_b3');

mean_sige1=mean(sige_b1');
std_sige1=std(sige_b1');
mean_sige2=mean(sige_b2');
std_sige2=std(sige_b2');
mean_sige3=mean(sige_b3');
std_sige3=std(sige_b3');

AD(1)=1;
AD=AD.*100;

AD=[1,10,100,1000];
xaxis=AD;

LW=2;MS=8;
%Fontsize
headF=13;
yaxisF=13;
fontF=11;

yminj=1;
ymaxj=2.5;

ymins=500;
ymaxs=2300;

yminse=500;
ymaxse=2300;

ymin=-150;
ymax=50;
XLABEL='Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)';
%XAXIS=[0.0507,0.1122,0.188,0.399,0.651,1.172];
XAXIS=[0.0507,0.1122,0.188,0.399,0.651];


XS=0.5;
XE=4.5;

xtix = {'100%','50%','25%','10%'};   % Your labels
xtixloc = [1 2 3 4];      % Your label locations


%plotting vars
axisF=23;
yaxisF=22;
 fontF=22;

hFig=figure;
h=errorbar(XAXIS,mean_J1,std_J1,'k','LineWidth',LW,'MarkerSize',MS);
%set (gca,'Xdir','reverse');
set([gca]             , ...
    'FontSize'   , fontF          );
ylabel('$J_{V}$','interpreter','latex','FontSize',yaxisF)
xlabel(XLABEL,'interpreter','latex','FontSize',yaxisF)
set(gca,'XTick', 0:0.2:0.8);
set(gca,'YTick', 1.05:0.05:1.2);
ylim([1.04 1.2]);
xlim([0 0.7]);

axis square
%set(gcf, 'PaperSize', [17 18])
%set(gcf, 'PaperPosition', [1.5 0.2 18 18]); %Position plot at left hand corner with width 5 and height 5.
output_plot_filename='j_Air_ball'

out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling-nmbe/figures/';

print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')


hFig=figure;
h=errorbar(XAXIS,mean_s_p1,std_s_p1,'k','LineWidth',LW,'MarkerSize',MS);
%set (gca,'Xdir','reverse');
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Pressure $(\mbox{Pa})$','interpreter','latex','FontSize',yaxisF)
  xlabel('Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)','interpreter','latex','FontSize',yaxisF)
set(gca,'XTick', 0:0.2:0.8);
 set(gca,'YTick', -150:50:0);
ylim([-190 0]);
xlim([0 0.7]);
%set(gcf, 'PaperSize', [16.1 15.5])
%set(gcf, 'PaperPosition', [-0.5 0.2 17 16]); %Position plot at left hand corner with width 5 and height 5.
axis square
output_plot_filename='p_Air_ball'
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')


hFig=figure;
h=errorbar(XAXIS,mean_sige1,std_sige1,'k','LineWidth',LW,'MarkerSize',MS);
%set (gca,'Xdir','reverse');
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Elastic stress $(\mbox{Pa})$','interpreter','latex','FontSize',yaxisF)
  xlabel('Pathway resistance ($\mbox{Pa}\,\mbox{mm}^{-3}\,\mbox{s}$)','interpreter','latex','FontSize',yaxisF)
set(gca,'XTick', 0:0.2:0.8);
 set(gca,'YTick', 550:50:750);
ylim([530 750]);
xlim([0 0.7]);
%set(gcf, 'PaperSize', [16.1 15.5])
%set(gcf, 'PaperPosition', [-0.5 0.2 17 16]); %Position plot at left hand corner with width 5 and height 5.
axis square

output_plot_filename='sige_Air_ball'
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')

