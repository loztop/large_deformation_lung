close all;
clear all;

Nu=289;
Np=512;

importfile('K.dat');
K = spconvert(data);
clear data;

importfile_rhs('r.dat');
r = sparse(data);
clear data;

Au=K(1:Nu,1:Nu);
Av=K(Nu+1:2*Nu,Nu+1:2*Nu);
Mx=K(2*Nu+Np+1:2*Nu+Np+Nu,2*Nu+Np+1:2*Nu+Np+Nu);
My=K(2*Nu+Np+Nu+1:2*Nu+Np+2*Nu,2*Nu+Np+Nu+1:2*Nu+Np+2*Nu);
J=K(2*Nu+1:2*Nu+Np,2*Nu+1:2*Nu+Np);
B_grad_u=K(1:Nu,2*Nu+1:2*Nu+Np);
B_grad_v=K(Nu+1:2*Nu,2*Nu+1:2*Nu+Np);
B_grad_x=K(2*Nu+Np+1:2*Nu+Np+Nu,2*Nu+1:2*Nu+Np);
B_grad_y=K(2*Nu+Np+Nu+1:2*Nu+Np+2*Nu,2*Nu+1:2*Nu+Np);
B_div_u=K(2*Nu+1:2*Nu+Np,1:Nu);
B_div_v=K(2*Nu+1:2*Nu+Np,Nu+1:2*Nu);
B_div_x=K(2*Nu+1:2*Nu+Np,2*Nu+Np+1:2*Nu+Np+Nu);
B_div_y=K(2*Nu+1:2*Nu+Np,2*Nu+Np+Nu+1:2*Nu+Np+2*Nu);



figure;
spy(K)



A=zeros(4*Nu+Np,4*Nu+Np);
A(1:Nu,1:Nu)=Au;
A(Nu+1:2*Nu,1+Nu:2*Nu)=Av;
A(2*Nu+1:3*Nu,1+2*Nu:3*Nu)=Mx;
A(3*Nu+1:4*Nu,1+3*Nu:4*Nu)=My;
A(4*Nu+1:4*Nu+Np,4*Nu+1:4*Nu+Np)=J;

A(1:Nu,4*Nu+1:4*Nu+Np)=B_grad_u;
A(Nu+1:2*Nu,4*Nu+1:4*Nu+Np)=B_grad_v;
A(2*Nu+1:3*Nu,4*Nu+1:4*Nu+Np)=B_grad_x;
A(3*Nu+1:4*Nu,4*Nu+1:4*Nu+Np)=B_grad_y;

A(4*Nu+1:4*Nu+Np,1:Nu)=B_div_u;
A(4*Nu+1:4*Nu+Np,Nu+1:2*Nu)=B_div_v;
A(4*Nu+1:4*Nu+Np,2*Nu +1:3*Nu)=B_div_x;
A(4*Nu+1:4*Nu+Np,3*Nu+1:4*Nu)=B_div_y;


fu=r(1:Nu,1);
fv=r(Nu+1:2*Nu,1);
fx=r(2*Nu+Np+1:2*Nu+Np+Nu,1);
fy=r(2*Nu+Np+Nu+1:2*Nu+Np+2*Nu,1);
fp=r(2*Nu+1:2*Nu+Np,1);

B=zeros(4*Nu+Np,1);
B(1:Nu,1)=fu;
B(Nu+1:2*Nu,1)=fv;
B(2*Nu+1:3*Nu,1)=fx;
B(3*Nu+1:4*Nu,1)=fy;
B(4*Nu+1:4*Nu+Np,1)=fp;

A=sparse(A)
B=sparse(B)

figure;
spy(A)

x = K\r;
sum(x)

x2 = A\B;
sum(x)


savefile = 'poro.mat';
save(savefile, 'Nu', 'Np','A','B','Au','Av','B_grad_u','B_grad_v','B_div_u','B_div_v','Mx','My','B_grad_x','B_grad_y','B_div_x','B_div_y','J','fu','fv','fy','fx','fy','fp')

