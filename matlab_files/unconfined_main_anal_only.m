clear all;


%Parameters for simulation (analytical solution)
params.T=100; %End Time
params.Nt=10;  %Number of time steps
params.dt=params.T/params.Nt;    %Size of time step

params.a=1;    %Radius of cylinder
params.v=0.15; %poisson ratio of elastic skeleton
params.E=10000; %Youngs modulus of elastic skeleton

params.lambda=(params.E*params.v)/((1+params.v)*(1-2*params.v));   %elastic coefficent
params.mu=params.E/(2*(1+params.v));   %elastic coefficent

%params.Hk=1;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)
params.Hk=params.lambda+2*params.mu;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)

params.k=0.000001;    %dynamic permeability
params.tg=1/(params.Hk*params.k/(params.a*params.a));  %characteristic time of diffusion
params.ez=0.05; %Amplitude of applied axial strain
params.Nodes=133; %Number of nodes

%create matrix to store results (NTx(number of idx x===1))
a=zeros(params.Nt,4);



%Calculate analytical solution
[b_y,b_x]=bessel(params);

%plot numerical solution
figure;
%plot analytical solution
plot(b_x,b_y,'k','LineWidth',3);
hold all



%axis([0 1.1 0 0.6])

% %Calculate Root-mean-square deviation of available co incideing time points
% rmsd=0;
% for i=1:params.Nt
% anal_y(i)=y(i*50);
% rmsd=rmsd+(anal_y(i)-num_y(i))^2;
% end
% rmsd=sqrt(rmsd/params.Nt)



