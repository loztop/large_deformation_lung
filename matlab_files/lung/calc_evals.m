% Given a real symmetric 3x3 matrix A, compute the eigenvalues
 A=magic(3)
 I=[1,0,0;0,1,0;0,0,1];
p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
if (p1 == 0) 
   % A is diagonal.
   eig1 = A(1,1)
   eig2 = A(2,2)
   eig3 = A(3,3)
else
   q = trace(A)/3
   p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
   p = sqrt(p2 / 6)
   B = (1 / p) * (A - q * I)       % I is the identity matrix
   r = det(B) / 2
 
   % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
   % but computation error can leave it slightly outside this range.
   if (r <= -1) 
      phi = pi / 3
   elseif (r >= 1)
      phi = 0
   else
      phi = acos(r) / 3
   end
 
   % the eigenvalues satisfy eig3 <= eig2 <= eig1
   eig1 = q + 2 * p * cos(phi)
   eig3 = q + 2 * p * cos(phi + (2*pi/3))
   eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
end