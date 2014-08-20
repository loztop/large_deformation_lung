clear all;
close all;

tec='.tec';

prefix='5036_Ce_';

%x_axis_str='Percentage of original Young''s modulus';
x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point/plot_data/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Crafeo_';

%home
base='/home/loztop/mount_point/mount_point/plot_data/';
out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Crafeo_';

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

C=[-60,-170,-90];

NP=size(P,1);

AD=[1,0.5,0.25,0.1,0.01];
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
 
        dist_C(j)= sqrt((C(1,1)-x(j)).^2 + (C(1,2)-y(j)).^2 + (C(1,3)-z(j)).^2);
  
    
end

%Pick out nodes for relevant bands of points

b1_c=0;
b2_c=0;
b3_c=0;

for j=1:Nodes
 
        if(dist_C(j)<22)
           b1_c=b1_c+1;
           b1_idx(b1_c)=j;
        end
  
        if(dist_C(j)>25 && dist_C(j)<27)
           b2_c=b2_c+1;
           b2_idx(b2_c)=j;
        end
        
        if(dist_C(j)>30 && dist_C(j)<40 )
           b3_c=b3_c+1;
           b3_idx(b3_c)=j;
        end
    
end


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

AD(1)=1;
AD=AD.*100;

AD=[1,10,100,1000];
xaxis=AD;

hFig=figure;

%Fontsize
headF=13;
yaxisF=13;
fontF=11;

yminj=0.98;
ymaxj=1.26;

ymins=1700;
ymaxs=4800;

yminse=150;
ymaxse=3100;

ymin=-2500;
ymax=200;

xtix = {'100%','50%','25%','10%','1%'};   % Your labels
xtixloc = [1 2 3 4 5];      % Your label locations

tight_subplot(4,3)
set(hFig, 'Position', [100 100 1000 1990])


set([gca]             , ...
    'FontSize'   , fontF           );

tight_subplot(4,3)
%%%plot Jacobian
h2=subplot(4,3,2)
h=boxplot(J_b2')
set(h(7,:),'Visible','off') 
hold all
 title('0-0.2cm from diseased area','interpreter','latex','FontSize',headF);
set([gca]             , ...
    'FontSize'   , fontF           );

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax2=get(h2,'Position');
ylim([yminj ymaxj])

 
subplot(4,3,3)
h=boxplot(J_b3')
set(h(7,:),'Visible','off') 
hold all
 title('0.5-1cm from diseased area','interpreter','latex','FontSize',headF);
set([gca]             , ...
    'FontSize'   , fontF           );


set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylim([yminj ymaxj])


h1=subplot(4,3,1)
h=boxplot(J_b1')
set(h(7,:),'Visible','off') 

hold all
 title('Inside diseased area','interpreter','latex','FontSize',headF);
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Jacobian','interpreter','latex','FontSize',yaxisF)

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax1=get(h1,'Position');
 set(h1,'Position',[ax1(1)-0.045 ax1(2) ax2(3) ax2(4)]);
ylim([yminj ymaxj])



%%%%%%%%%%%%%%%%%%%%%%%%
h5=subplot(4,3,5)
h=boxplot(s_p_b2')
set(h(7,:),'Visible','off') 
hold all
% title('0-0.2cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax5=get(h5,'Position');
ylim([ymin ymax])

 
subplot(4,3,6)
h=boxplot(s_p_b3')
set(h(7,:),'Visible','off') 
hold all
% title('0.5-1cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );


set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylim([ymin ymax])


h4=subplot(4,3,4)
h=boxplot(s_p_b1')
set(h(7,:),'Visible','off') 
hold all
% title('Inside diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Pressure (Pa)','interpreter','latex','FontSize',yaxisF)

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax4=get(h4,'Position');
 set(h4,'Position',[ax4(1)-0.063 ax4(2) ax5(3) ax5(4)]);
ylim([ymin ymax])


%%%%%%%%%%%%%%%%%%%%%%%% STRESS
h8=subplot(4,3,8)
h=boxplot(sig_b2')
set(h(7,:),'Visible','off') 
hold all
% title('0-0.2cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax8=get(h8,'Position');
ylim([ymins ymaxs])

 
subplot(4,3,9)
h=boxplot(sig_b3')
set(h(7,:),'Visible','off') 
hold all
% title('0.5-1cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );


set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylim([ymins ymaxs])


h7=subplot(4,3,7)
h=boxplot(sig_b1')
set(h(7,:),'Visible','off') 
hold all
% title('Inside diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Total stress (Pa)','interpreter','latex','FontSize',yaxisF)

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax7=get(h7,'Position');
 set(h7,'Position',[ax7(1)-0.055 ax7(2) ax8(3) ax8(4)]);
ylim([ymins ymaxs])
 

%%%%%%%%%%%%%%%%%%%%%%%% ELASTIC STRESS
h8=subplot(4,3,11)
h=boxplot(sige_b2')
set(h(7,:),'Visible','off') 
hold all
% title('0-0.2cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax8=get(h8,'Position');
ylim([yminse ymaxse])

 
subplot(4,3,12)
h=boxplot(sige_b3')
set(h(7,:),'Visible','off') 
hold all
% title('0.5-1cm from diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );


set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylim([yminse ymaxse])


h7=subplot(4,3,10)
h=boxplot(sige_b1')
set(h(7,:),'Visible','off') 
hold all
% title('Inside diseased area','interpreter','latex','FontSize',14);
set([gca]             , ...
    'FontSize'   , fontF           );
  ylabel('Elastic stress (Pa)','interpreter','latex','FontSize',yaxisF)

set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ax7=get(h7,'Position');
 set(h7,'Position',[ax7(1)-0.055 ax7(2) ax8(3) ax8(4)]);
ylim([yminse ymaxse])

 set(gcf, 'PaperPosition', [-1 -2.2 22 27]); %Position plot at left hand corner with width 5 and height 5.
 set(gcf, 'PaperSize', [20 24]);

 
output_plot_filename='5036'
print(hFig,strcat(out_base,output_plot_filename),'-dpdf','-r0')
 

