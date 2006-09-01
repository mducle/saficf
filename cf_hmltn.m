function M = cf_hmltn(J,B2,B4,B6)
% cf_hmltn - Calculates the crystal field Hamiltonian under Steven's normalisation. 
%
% Syntax:  Hcf = cf_hmltn(J,B2,B4,B6)
%     OR:  Hcf = cf_hmltn(J,{B2 B4 B6})  [where the parameters are in a cell array]
%
% Inputs:  J   - The total angular momentum quantum number 
%          B2  - A 5-component vector containing the empirical crystal field 
%                parameters of the magnetic ion of interest of rank 2
%          B4  - is a 9-component vector with the parameters of rank 4
%          B6  - is a 13-component vector with the parameters of rank 6
%                The format for each vector is from B_{-(rank)} to B_{+(rank)}
%
% Outputs: Hcf - A (2J+1)x(2J+1) matrix representing the Crystal Field Hamiltonian
%
% This function is based on the program ENGYLVL.BAS by Newman and Ng
% From their book: Crystal Field Handbook (CUP 2000).

% By Duc Le 2005 - duc.le@ucl.ac.uk
% mdl - updated 060730: For saficf, improved some clunky contructions. 
%     - updated 060815: Added cell-array input for norm_cfpars use.
%     - updated 060820: Backported faster routines from norm_cfhmltn.m

% This file is part of the SAfiCF package. 
% Licenced under the GNU GPL v2 or later. 

% Checks that we have the right arguments.
if nargin == 2
  if ~iscell(B2)
    error(['If using only two arguments, B must be a cell array. Type "help ' ...
           mfilename '".']);
  else
    B4 = B2{2}; B6 = B2{3}; B2 = B2{1};
    if length(B2) ~= 5 || length(B4) ~= 9 || length(B6) ~=13
      error(['Invalid crystal field parameters. Type "help ' mfilename '".'])
    end
  end
elseif nargin == 4
  if ~isvector(B2) || length(B2) ~= 5
    error(['Invalid crystal field parameters B2. Type "help ' mfilename '".'])
  elseif ~isvector(B4) || length(B4) ~= 9
    error(['Invalid crystal field parameters B4. Type "help ' mfilename '".'])
  elseif ~isvector(B6) || length(B6) ~= 13
    error(['Invalid crystal field parameters B6. Type "help ' mfilename '".'])
  end
else
  error(['This function requires either four or two arguments. Type "help ' ...
         mfilename '".'])
end
% For backwards compatibility. Old version of cf_hmltn.m required all the angular 
% momentum quantum numbers for calculations of gJ and induced moments.
if isvector(J) && length(J) == 3
  L = J(1); S = J(2); J = J(3);
end
if ~isscalar(J) || mod((J*2),1)
  error(['Total angular momentum J must be a scalar and integer or half integer.' ...
         ' Type "help ' mfilename '".'])
end

% The crystal field potential operator is given by:
%
%         ---   k    [  k        q  k    ]     ---  k  k     ---   k [  k       q  k  ]
% V   = i >    B     | O   - (-1)  O     |  +  >   B  O   +  >    B  | O  + (-1)  O   |
%  cf     ---   -|q| [  |q|         -|q| ]     ---  0  0     ---   q [  q          -q ]
%        k,q<0                                  k           k,q>0  
%
% where O_k^q are the Stevens operators of rank k.
% Reference: C. Rudowicz, J. Phys. C: Solid State Phys., vol 18, pp1415-1430 (1985).
%
% The energy matrix elements <LSJM_j|V_CF|L'SJ'M_j'> are given by summing
% over the matrix elements of the Stevens operator given by:
%
%          k                 J-M_j                           k
% <LSJM | O  |L'SJ'M'> = (-1)      ( J  k  J' )   *   (LSJ||O ||L'SJ')
%      j   q        j              (-Mj q  Mj')
%
%   where this is the 3-j symbol---^  and reduced matrix element--^
%
% The reduced matrix elements for the Stevens operators are given by:
%                   _______________
%      k        -k | (2J + k + 1)!
% <j||O ||j> = 2   | -------------
%                 \|   (2J - k)!
%
% Reference: D.Smith and J.H.M. Thornley, Proc. Phys. Soc., 1966, vol 89, pp779.

if 2*J-2>0
  RM2 = (1/4) * sqrt( factorial(2*J+2+1) / factorial(2*J-2) );
else 
  RM2 = 0;
end
if 2*J-4>0
  RM4 = (1/16) * sqrt( factorial(2*J+4+1) / factorial(2*J-4) );
else
  RM4 = 0;
end
if 2*J-6>0
  RM6 = (1/64) * sqrt( factorial(2*J+6+1) / factorial(2*J-6) );
else
  RM6 = 0;
end

% Calculates the multiplicity
dimj = 2*J+1;

% Defines the matrix
M = zeros(dimj);

% Hamiltonian must be hermitian, so M = M' (i.e. Mij = conj(Mji) ). So we only need to
% compute values of the upper (or lower) triangle.

% Each order q operator has only nonzero matrix elements along a diagonal of the matrix
% Eg. Order 0 terms has it along the diagonal. Order 1 terms has it one element left,
% and so on. We can use this to calculate only for the orders with nonzero parameters.

% NB. the matrix elements calculated here using the threej symbols are actually matrix 
% elements of the Wybourne normalised operator equivalent to spherical harmonics rather
% than the Steven's operator equivalent to tesseral harmonics. However both sets will 
% give the same eigenvector and eigenvalues when diagonalised. Nevertheless, in order to
% get numerical agreements, we have to multiply each rank-k, component-q term by the 
% ratios of the Wybourne to the Stevens parameters, which is the factor at the end of
% each term in the sums. See the Rudowicz paper for clarification.

% Zero component terms: B20, B40, B60
if sum(find(abs(find(B2)-3)==0)) | sum(find(abs(find(B4)-5)==0)) | sum(find(abs(find(B6)-7)==0))
  for i = 1:dimj
    Mj = i-1 - J;
    % Rank 2
    M(i,i) = M(i,i) + (-1)^(J-Mj) * threej([J 2 J; -Mj  0 Mj]) * RM2 * B2(3) * 2;    
    % Rank 4
    M(i,i) = M(i,i) + (-1)^(J-Mj) * threej([J 4 J; -Mj  0 Mj]) * RM4 * B4(5) * 8; 
    % Rank 6
    M(i,i) = M(i,i) + (-1)^(J-Mj) * threej([J 6 J; -Mj  0 Mj]) * RM6 * B6(7) * 16;
  end
end % if nonzero zeroth component.

% One component terms: B21, B2-1, B41, B4-1, B61, B6-1
if sum(find(abs(find(B2)-3)==1)) | sum(find(abs(find(B4)-5)==1)) | sum(find(abs(find(B6)-7)==1))
  for i = 1:dimj-1 
    j = i + 1;
    Mj = i-1 - J;
    Mjp = i - J;
    % Calculates the threej values here to save computation time later.
    tj21 = threej([J 2 J; -Mj  1 Mjp]);     tj2m1 = threej([J 2 J; -Mj -1 Mjp]);
    tj41 = threej([J 4 J; -Mj  1 Mjp]);     tj4m1 = threej([J 4 J; -Mj -1 Mjp]);
    tj61 = threej([J 6 J; -Mj  1 Mjp]);     tj6m1 = threej([J 6 J; -Mj -1 Mjp]);
    % Rank 2
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj21  * RM2 * B2(2) *-sqrt(-6);
    M(i,j) = M(i,j) - (-1)^(J-Mj-1) * tj2m1 * RM2 * B2(2) *-sqrt(-6);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj21  * RM2 * B2(4) *-sqrt(6);
    M(i,j) = M(i,j) + (-1)^(J-Mj+1) * tj2m1 * RM2 * B2(4) *-sqrt(6);
    % Rank 4
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj41  * RM4 * B4(4) *-2/sqrt(-5);
    M(i,j) = M(i,j) - (-1)^(J-Mj-1) * tj4m1 * RM4 * B4(4) *-2/sqrt(-5);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj41  * RM4 * B4(6) *-2/sqrt(5);
    M(i,j) = M(i,j) + (-1)^(J-Mj+1) * tj4m1 * RM4 * B4(6) *-2/sqrt(5);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj61  * RM6 * B6(6) *-8/sqrt(-42);
    M(i,j) = M(i,j) - (-1)^(J-Mj-1) * tj6m1 * RM6 * B6(6) *-8/sqrt(-42);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj61  * RM6 * B6(8) *-8/sqrt(42);
    M(i,j) = M(i,j) + (-1)^(J-Mj+1) * tj6m1 * RM6 * B6(8) *-8/sqrt(42);
  end
end % if nonzero one component.

% Two component terms: B22, B2-2, B42, B4-2, B62, B6-2
if sum(find(abs(find(B2)-3)==2)) | sum(find(abs(find(B4)-5)==2)) | sum(find(abs(find(B6)-7)==2))
  for i = 1:dimj-2 
    j = i + 2;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj22 = threej([J 2 J; -Mj  2 Mjp]);     tj2m2 = threej([J 2 J; -Mj -2 Mjp]);
    tj42 = threej([J 4 J; -Mj  2 Mjp]);     tj4m2 = threej([J 4 J; -Mj -2 Mjp]);
    tj62 = threej([J 6 J; -Mj  2 Mjp]);     tj6m2 = threej([J 6 J; -Mj -2 Mjp]);
    % Rank 2
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj22  * RM2 * B2(1) * 2/sqrt(-6);
    M(i,j) = M(i,j) - (-1)^(J-Mj-2) * tj2m2 * RM2 * B2(1) * 2/sqrt(-6);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj22  * RM2 * B2(5) * 2/sqrt(6);
    M(i,j) = M(i,j) + (-1)^(J-Mj+2) * tj2m2 * RM2 * B2(5) * 2/sqrt(6);
    % Rank 4
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj42  * RM4 * B4(3) * 4/sqrt(-10);
    M(i,j) = M(i,j) - (-1)^(J-Mj-2) * tj4m2 * RM4 * B4(3) * 4/sqrt(-10);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj42  * RM4 * B4(7) * 4/sqrt(10);
    M(i,j) = M(i,j) + (-1)^(J-Mj+2) * tj4m2 * RM4 * B4(7) * 4/sqrt(10);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj62  * RM6 * B6(5) * 16/sqrt(-105);
    M(i,j) = M(i,j) - (-1)^(J-Mj-2) * tj6m2 * RM6 * B6(5) * 16/sqrt(-105);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj62  * RM6 * B6(9) * 16/sqrt(105);
    M(i,j) = M(i,j) + (-1)^(J-Mj+2) * tj6m2 * RM6 * B6(9) * 16/sqrt(105);
  end
end % if nonzero two component

% Three component terms: B43, B4-3, B63, B6-3
if sum(find(abs(find(B4)-5)==3)) | sum(find(abs(find(B6)-7)==3))
  for i = 1:dimj-3 
    j = i + 3;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj43 = threej([J 4 J; -Mj  3 Mjp]);     tj4m3 = threej([J 4 J; -Mj -3 Mjp]);
    tj63 = threej([J 6 J; -Mj  3 Mjp]);     tj6m3 = threej([J 6 J; -Mj -3 Mjp]);
    % Rank 4
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj43  * RM4 * B4(2) *-2/sqrt(-35);
    M(i,j) = M(i,j) - (-1)^(J-Mj-3) * tj4m3 * RM4 * B4(2) *-2/sqrt(-35);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj43  * RM4 * B4(8) *-2/sqrt(35);
    M(i,j) = M(i,j) + (-1)^(J-Mj+3) * tj4m3 * RM4 * B4(8) *-2/sqrt(35);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj63  * RM6 * B6(4) *-8/sqrt(-105);
    M(i,j) = M(i,j) - (-1)^(J-Mj-3) * tj6m3 * RM6 * B6(4) *-8/sqrt(-105);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj63  * RM6 * B6(10) *-8/sqrt(105);
    M(i,j) = M(i,j) + (-1)^(J-Mj+3) * tj6m3 * RM6 * B6(10) *-8/sqrt(105);
  end
end % if nonzero three component

% Four component terms: B44, B4-4, B64, B6-4
if sum(find(abs(find(B4)-5)==4)) | sum(find(abs(find(B6)-7)==4))
  for i = 1:dimj-4 
    j = i + 4;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj44 = threej([J 4 J; -Mj  4 Mjp]);     tj4m4 = threej([J 4 J; -Mj -4 Mjp]);
    tj64 = threej([J 6 J; -Mj  4 Mjp]);     tj6m4 = threej([J 6 J; -Mj -4 Mjp]);
    % Rank 4
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj44  * RM4 * B4(1) * 8/sqrt(-70);
    M(i,j) = M(i,j) - (-1)^(J-Mj-4) * tj4m4 * RM4 * B4(1) * 8/sqrt(-70);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj44  * RM4 * B4(9) * 8/sqrt(70);
    M(i,j) = M(i,j) + (-1)^(J-Mj+4) * tj4m4 * RM4 * B4(9) * 8/sqrt(70);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj64  * RM6 * B6(3) * 16/3/sqrt(-14);
    M(i,j) = M(i,j) - (-1)^(J-Mj-4) * tj6m4 * RM6 * B6(3) * 16/3/sqrt(-14);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj64  * RM6 * B6(11) * 16/3/sqrt(14);
    M(i,j) = M(i,j) + (-1)^(J-Mj+4) * tj6m4 * RM6 * B6(11) * 16/3/sqrt(14);
  end
end % if nonzero four component.

% Five component terms: B65, B6-5
if sum(find(abs(find(B6)-7)==5))
  for i = 1:dimj-5 
    j = i + 5;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj65 = threej([J 6 J; -Mj  5 Mjp]);     tj6m5 = threej([J 6 J; -Mj -5 Mjp]);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj65  * RM6 * B6(2) *-8/3/sqrt(-77);
    M(i,j) = M(i,j) - (-1)^(J-Mj-5) * tj6m5 * RM6 * B6(2) *-8/3/sqrt(-77);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj65  * RM6 * B6(12) *-8/3/sqrt(77);
    M(i,j) = M(i,j) + (-1)^(J-Mj+5) * tj6m5 * RM6 * B6(12) *-8/3/sqrt(77);
  end
end % if nonzero five component.

% Six component terms: B66, B6-6
if sum(find(abs(find(B6)-7)==6))
  for i = 1:dimj-6
    j = i + 6;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj66 = threej([J 6 J; -Mj  6 Mjp]);     tj6m6 = threej([J 6 J; -Mj -6 Mjp]);
    % Rank 6
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj66  * RM6 * B6(1) * 16/sqrt(-231);
    M(i,j) = M(i,j) - (-1)^(J-Mj-6) * tj6m6 * RM6 * B6(1) * 16/sqrt(-231);
    M(i,j) = M(i,j) + (-1)^(J-Mj)   * tj66  * RM6 * B6(13) * 16/sqrt(231);
    M(i,j) = M(i,j) + (-1)^(J-Mj+6) * tj6m6 * RM6 * B6(13) * 16/sqrt(231);
  end
end % if nonzero six component

% Calculates the lower triangle
LT = M'; LT(1:dimj+1:dimj^2) = 0;

% Outputs the Hamiltonian.
M = M + LT;
