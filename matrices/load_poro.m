%close all;
%clear all;
load poro.mat;
whos

%Nu Number of p1 dofs
%Np Number of p0 dofs

%   A                                               B

%   Au      0       0       0          B_grad_u     fu
%   0       Av      0       0          B_grad_v     fv
%   0       0       Mx      0          B_grad_x     fx
%   0       0       0       0          B_grad_y     fy
%   B_div_u B_div_v B_div_x B_div_y    J            fb

figure;
spy(A);

x=A\B;