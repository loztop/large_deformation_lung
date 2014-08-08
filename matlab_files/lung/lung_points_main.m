clear all;
%close all;

tec='.tec';
G=0.435
l=0.186
E=G*(3*l+2*G)/(l+G)
v=l/(2*(l+G))

%numerical solution

base='/auto/users/lorenzb/mount_point3/whole_lung_large_deformation/data/';



Nodes=246 %for 246 mesh

%Tracking points
P1=[75,110,160];
P2=[105,110,160];
P3=[135,110,160];

P=[75,110,160;
    105,110,160;
    135,110,160];

NP=size(P,1);

AD=[0,0.5,0.25,0.01];
TD=[0,0,0,0];

NT=1;

for i=1:length(AD)
i
%Load time step data    
fnameC=strcat(base,'C_',num2str(AD(i)),'_',num2str(TD(i)),'__',num2str(NT),tec);
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
near_P_idx;

  
%The displacements
%s_u=data(:,4);
%s_v=data(:,5);
%s_w=data(:,6);


%Pressure
s_p(i,:)=data(near_P_idx(i,:),7)';
%Jacobian
J(i,:)=data(near_P_idx(i,:),31)';
%'p1' - the toatal stress
sig(i,:)=data(near_P_idx(i,:),32)';

end

AD(1)=1;
figure;
plot(AD,s_p(:,1),'-x');
hold all;
plot(AD,s_p(:,2),'-o');
hold all;
plot(AD,s_p(:,3),'-.');
hold all;
title('Aveolar pressure');

figure;
plot(AD,J(:,1),'-x');
hold all;
plot(AD,J(:,2),'-o');
hold all;
plot(AD,J(:,3),'-.');
hold all;
title('Jacobian');

figure;
plot(AD,sig(:,1),'-x');
hold all;
plot(AD,sig(:,2),'-o');
hold all;
plot(AD,sig(:,3),'-.');
hold all;
title('sigma');



%plot numerical solution
%plot(num_x,num_y,'rx','MarkerSize',14,'LineWidth',2);
%hold all
%save matfiles/unconfined_non_p001_726_250.mat


