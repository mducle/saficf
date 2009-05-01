function Ji = mag_op_j(J)
% mag_op_j - returns the matrix representations for operators Jx, Jy, Jz
%
% Syntax:  Ji = mag_op_j(J)
%
% Inputs:  J is a scalar being the total angular momentum quantum number
%
% Outputs: Ji is a (2J+1)x(2J+1)x5 matrix array, with Ji(:,:,i) where i=x,y,z,+,-

if size(J,2) ~= 1
  error('Require J scalar');
end

% The dipole magnetic operators Jx, Jy, Jz are equivalent to rank 1 (k=1)
% Stevens operators, and are given by:
%
%          1    -1                  1     -1              0 
% J  = -( O  - O   )      J  =  i( O  +  O   )      J  = O 
%  x       1    1          y        1     1          z    1
% 
% Reference: K.W.H. Stevens, Proc. Phys. Soc. A, 1952,  vol 65, pp209.

% The matrix elements of the Stevens operator given by:
%                                                         ____________
%          k                 J-M_j                     -k | (2J+k+1)!
% <LSJM | O  |L'SJ'M'> = (-1)      ( J  k  J' )   *   2   | ---------
%      j   q        j              (-Mj q  Mj')          \|  (2J-k)!
%
%   where this is the 3-j symbol---^  and reduced matrix element--^
%
% Reference: D.Smith and J.H.M. Thornley, Proc. Phys. Soc., 1966, vol 89, pp779.

% Works out the reduced matrix elements
%if 2*J-1>0
  RM1 = (1/2) * sqrt( factorial(2*J+1+1) / factorial(2*J-1) );
%else 
%  RM1 = 1;
%end

ind_i = 0;
for Mj = -J:J
  ind_i = ind_i + 1;
  ind_j = 0;
  for Mjp = -J:J
    ind_j = ind_j + 1;

% Rank 1
    Om(ind_i,ind_j) = (-1)^(J-Mj) * threej([J 1 J; -Mj -1 Mjp]) * RM1 / sqrt(2);
    O0(ind_i,ind_j) = (-1)^(J-Mj) * threej([J 1 J; -Mj  0 Mjp]) * RM1 * 1;
    Op(ind_i,ind_j) = (-1)^(J-Mj) * threej([J 1 J; -Mj  1 Mjp]) * RM1 / sqrt(2);

  end
end

Ji(:,:,1) =  -( Op - Om );   % Jx
Ji(:,:,2) = i*( Op + Om );   % Jy
Ji(:,:,3) =     O0;          % Jz
Ji(:,:,4) = -Op;
Ji(:,:,5) = Om;

