clear all;
close all;

tec='.tec';

prefix='1447_Cr_';
AD=[1,0.6,0.5,0.4];
TD=[0,0,0,0];

%prefix='5036_Wr_';
%AD=[0];
%TD=[0.1];

x_axis_str='Percentage of original Young''s modulus';
%x_axis_str='Percentage of original airway constriction';

%numerical solution

%work
base='/auto/users/lorenzb/mount_point/plot_data/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Wrafeo_';

%home
base='/home/loztop/mount_point/mount_point/plot_data/';
out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/Wrafeo_';

Nodes=1447; %for 246 mesh

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

for i=length(AD):length(AD)
i
%Load time step data    
fnameC=strcat(base,prefix,num2str(AD(i)),'_',num2str(TD(i)),'__',num2str(NT),tec)
importfile(fnameC);
header_info=textdata{3};
TEMP=data;
clear data;
dataD=TEMP(1:Nodes,:);


fnameH=strcat(base,prefix,num2str(1),'_',num2str(0),'__',num2str(NT),tec)
importfile(fnameH);
header_info=textdata{3};
TEMP=data;
clear data;
dataH=TEMP(1:Nodes,:);


%This might change depending on what order vraibles are stored

t=size(dataD,1);

%reference positions
x=dataD(:,20)
y=dataD(:,21);
z=dataD(:,22);

%loop thorough the nodes to find the closest one to the three points of
%interest

near_dist_P=ones(1,4).*999;
near_P_idx(i,:)=zeros(1,4);
for j=1:Nodes
 
        dist_C(j)= sqrt((C(1,1)-x(j)).^2 + (C(1,2)-y(j)).^2 + (C(1,3)-z(j)).^2);
  
    
end

%Pick out nodes for relevant bands of points

DIST=[16:1:39];
b_c=zeros(length(DIST)-1,Nodes);

for j=1:Nodes
 
    for d=1:length(DIST)-1
        if(dist_C(j)>DIST(d) && dist_C(j)<DIST(d+1))
           b_c(d)=b_c(d)+1;
           b_idx(d,b_c(d))=j;
        end
    end
  j
end

%%Jacobian

i
end

 for d=1:length(DIST)-1
     mJ(d)=mean(dataD(b_idx(d),34))';
     sJ(d)=std(dataD(b_idx(d),34))';

     mP(d)=mean(dataD(b_idx(d),7))';
     sP(d)=std(dataD(b_idx(d),7))';
     
     mSt(d)=mean(dataD(b_idx(d),33))';
     sSt(d)=std(dataD(b_idx(d),33))';
     
     mSe(d)=mean(dataD(b_idx(d),32))';
     sSe(d)=std(dataD(b_idx(d),32))';
     
     %Healthy data
     mJh(d)=mean(dataH(b_idx(d),34))';
     sJh(d)=std(dataH(b_idx(d),34))';
     
     mPh(d)=mean(dataH(b_idx(d),7))';
     sPh(d)=std(dataH(b_idx(d),7))';
     
     mSeh(d)=mean(dataH(b_idx(d),32))';
     sSeh(d)=std(dataH(b_idx(d),32))';
     
 end
 
 
 mean_sigz_A(i)=median(SiA);
std_sigz_A(i)=std(SiA);
sigz_AQ1(i) = median(SiA(find(SiA<median(SiA))));
sigz_AQ3(i) = median(SiA(find(SiA>median(SiA))));   
 
 
 figure;
 plot(DIST(1:end-1),mJ);
 hold all;
 plot(DIST(1:end-1),mJh,'r');

 figure;
 plot(DIST(1:end-1),mP);
 hold all;
 plot(DIST(1:end-1),mPh,'r');
 
 figure;
 plot(DIST(1:end-1),mSe);
 hold all;
 plot(DIST(1:end-1),mSeh,'r');