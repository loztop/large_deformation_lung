clear all;
close all;

tec='.tec';

base='/auto/users/lorenzb/mount_point/data/';



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


%Load time step data    
fnameC=strcat(base,'Health80NT4T_20',tec);
fnameC=strcat(base,'Const1339_0_6NT2T_6',tec);
fnameC=strcat(base,'N48_2881_0.01E_fix_8T_80NT_10',tec);

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


%Pressure
s_p=data(:,7)';
%Jacobian
J=data(:,31)';
%'p1' - the toatal stress
sig=data(:,32)';


mean(J)
mean(sig)
mean(s_p)




