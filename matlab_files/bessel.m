function [y,x] = bessel(params)

format long

%Load parameters
v=params.v;
a=params.a;
k=params.k;
ez=params.ez;
Hk=params.Hk;
T=params.T;
tg=params.tg;

[alpha] = find_roots(v); %Roots of bessel function

Nt=20000; %Number of time steps for analytical solution
%Nt=params.Nt; %Number of time steps for analytical solution

u=zeros(1,Nt);  %Solution vector of analytical radial displacement
dt=T/Nt;
t=0;


for j=1:Nt
	summ=0;
    for i=1:length(alpha)	  
        top=exp(-(alpha(i)^2)*(1/tg)*t);
        bot=(alpha(i)^2)*((1-v)^2)-(1-2*v);
        summ=summ+top/bot;
    end
	u(j)=a*ez*(v + (1-2*v)*(1-v)*summ);
    t=t+dt;
end


%Output x axis vector and solution (scaled)
y=u/(a*ez);
x=dt*(1:1:Nt)./tg;

%Plot figure of analytical solution
%figure;
%plot(x,y);
