clear all;
close all;

tec='.tec';

%work
%base='/auto/users/lorenzb/mount_point3/whole_lung_large_deformation/data/';

%home
base='/home/loztop/mount_point2/mount_point3/whole_lung_large_deformation/data/';


Nodes=246 %for 246 mesh

%Tracking points
 
P=[75,110,160;
    75,130,152;
    75,152,144];

NP=size(P,1);

TD=[0,0.5,0.25,0.1,0.01];
AD=[0,0,0,0,0];

NT=2;

for i=1:length(AD)
i
%Load time step data    
fnameC=strcat(base,'W_',num2str(AD(i)),'_',num2str(TD(i)),'__',num2str(NT),tec);
importfile(fnameC);

fnameC

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

near_dist_P=[999,999,999];
near_P_idx(i,:)=[0,0,0];
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

TD(1)=1;
figure;
plot(TD,s_p(:,1),'-x');
hold all;
plot(TD,s_p(:,2),'-o');
hold all;
plot(TD,s_p(:,3),'-.');
hold all;
title('Aveolar pressure');

figure;
plot(TD,J(:,1),'-x');
hold all;
plot(TD,J(:,2),'-o');
hold all;
plot(TD,J(:,3),'-.');
hold all;
title('Jacobian');

figure;
plot(TD,sig(:,1),'-x');
hold all;
plot(TD,sig(:,2),'-o');
hold all;
plot(TD,sig(:,3),'-.');
hold all;
title('sigma');



%plot numerical solution
%plot(num_x,num_y,'rx','MarkerSize',14,'LineWidth',2);
%hold all
%save matfiles/unconfined_non_p001_726_250.mat


