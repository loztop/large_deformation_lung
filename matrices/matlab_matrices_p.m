close all;
%clear all;


importfile('Kp_tn2000.dat');
Kz = spconvert(data);
clear data;

condest(Kz)


importfile('Kp_pn2000.dat');
Kx = spconvert(data);
clear data;

condest(Kx)

importfile('Kp_pnsym2000.dat');
Kpnsym = spconvert(data);
clear data;

condest(Kpnsym)


importfile_rhs('rp.dat');
rp = sparse(data);
clear data;


figure;
spy(Kz)


figure;
spy(Kx)

figure;
spy(Kpnsym)

%figure;
%spy(Kpt)



