clear all;
close all;


tec='.tec';


prefix='plot_data/5036_Cr_1_0__4';
Nodes=5036; 
FRC=1.477; %1447
REF=0.907; %1447
 
prefix_res='data/5036_H_air_res_0';
prefix_res='data/5036_C_0.35_air_res_0';


NT=[0:1:20];


%work
base='/auto/users/lorenzb/mount_point/';
out_base='/users/lorenzb/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';

%home
base='/home/loztop/mount_point/mount_point/';
%out_base='/home/loztop/Dropbox/Dphil/poroelasticity_papers/coupling_paper/figures/lung_sims/';



 i=1;


%get airway_tree resistance
fnameC=strcat(base,prefix_res,tec)
importfile(fnameC);
header_info=textdata{3};
TEMP=data;
clear data;
data_res=TEMP(1:Nodes,:);
R=data_res(:,31)';

%reference positions
x=data_res(:,20);
y=data_res(:,21);
z=data_res(:,22);

%loop thorough the nodes to find the closest one to the three points of
%interest
C=[-57,-144,-90];

near_dist_P=ones(1,4).*999;
near_P_idx(i,:)=zeros(1,4);
for j=1:Nodes
 
        dist_C(j)= sqrt((C(1,1)-x(j)).^2 + (C(1,2)-y(j)).^2 + (C(1,3)-z(j)).^2);
  
    
end

b1_c=0;
for j=1:Nodes
 
        if(dist_C(j)<19)
           b1_c=b1_c+1;
           b_idx(b1_c)=j;
        end
    
end

mean(R(b_idx))

