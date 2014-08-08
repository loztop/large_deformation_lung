clear all;
%close all;

tec='.tec';
G=0.435
l=0.186
E=G*(3*l+2*G)/(l+G)
v=l/(2*(l+G))

%numerical solution

base='/auto/users/lorenzb/mount_point/examples/poroelasticity/nonlinear_poroelasticity/data/cylinder_sym138_20NT_10T_p001_STAB_';
base='/auto/users/lorenzb/mount_point/examples/poroelasticity/nonlinear_poroelasticity/data/cylinder_sym728_20NT_10T_p001_STAB_';
base='/auto/users/lorenzb/mount_point/examples/poroelasticity/nonlinear_poroelasticity/data/cylinder_sym5567_100NT_10T_p001_STAB_';

base='/auto/users/lorenzb/mount_point/examples/poroelasticity/nonlinear_poroelasticity/data/cylinder_sym728_250NT_1000T_p001_STAB_';

%base='/auto/users/lorenzb/mount_point/examples/poroelasticity/nonlinear_poroelasticity/data/cylinder_sym728_100NT_1000T_p001_STAB_';


%Parameters for simulation (analytical solution)
params.T=1000; %End Time
params.Nt=250;  %Number of time steps
params.Nodes=726; %Number of nodes  5565, 726

params.dt=params.T/params.Nt;    %Size of time step

params.a=1;    %Radius of cylinder
params.v=0.15; %poisson ratio of elastic skeleton
params.E=1; %Youngs modulus of elastic skeleton

params.lambda=(params.E*params.v)/((1+params.v)*(1-2*params.v));   %elastic coefficent
params.mu=params.E/(2*(1+params.v));   %elastic coefficent

%params.Hk=1;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)
params.Hk=params.lambda+2*params.mu;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)

params.k=0.001;    %dynamic permeability
params.tg=1/(params.Hk*params.k/(params.a*params.a));  %characteristic time of diffusion
params.ez=0.01; %Amplitude of applied axial strain

%create matrix to store results (NTx(number of idx x===1))
a=zeros(params.Nt,4);


for i=1:params.Nt
i
%Load time step data    
fname=strcat(base,num2str(i),tec);
importfile(fname);

header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:params.Nodes,:);

%This might change depending on what order vraibles are stored

t=size(data,1);

%The reference positions
x=data(:,18);
y=data(:,19);
z=data(:,20);

%The displacements
s_u=data(:,4);
s_v=data(:,5);
s_w=data(:,6);

%Find index of nodes where x=1, y=1, z=1
idx_x_1=find(x==1);
idx_y_1=find(y==1);
idx_z_1=find(z==1);

%get the radial (x direction) displacemnet at points x=1
s_u_x_1=s_u(idx_x_1);

%Store these displacement for this time in a vector
%a(i,:)=s_u_x_1;

%Store the max and min values of radial displacement
%num_u(i)=max(s_u);
%num_u_min(i)=min(s_u);

num_u_min(i)=min(s_v);

num_u(i)=(max(s_v)-min(s_v))/2;


end
num_u

%Creat numeical x-axis and y solution vector
num_x=params.dt.*(1:5:params.Nt)./params.tg;
num_y=1.*(num_u(1:5:end)-1)./(1*params.a*params.ez);

%Calculate analytical solution
[b_y,b_x]=bessel(params);

figure;
%plot analytical solution
plot(b_x,b_y,'k','LineWidth',3);
hold all

%plot numerical solution
plot(num_x,num_y,'rx','MarkerSize',14,'LineWidth',2);
hold all


save matfiles/unconfined_non_p001_726_250.mat
%clear
%load unconfined.mat

%axis([0 1.1 0 0.6])

% %Calculate Root-mean-square deviation of available co incideing time points
% rmsd=0;
% for i=1:params.Nt
% anal_y(i)=y(i*50);
% rmsd=rmsd+(anal_y(i)-num_y(i))^2;
% end
% rmsd=sqrt(rmsd/params.Nt)



