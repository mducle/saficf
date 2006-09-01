function M = norm_cfhmltn(J,V2,V4,V6)
% cf_hmltn - Calculates the crystal field Hamiltonian under Fabi normalisation. 
%
% Syntax:  Hcf = cf_hmltn(J,V2,V4,V6)
%     OR:  Hcf = cf_hmltn(J,{V2 V4 V6})  [where the parameters are in a cell array]
%
% Inputs:  J   - The total angular momentum quantum number 
%          V2  - A 5-component vector containing the empirical crystal field 
%                parameters of the magnetic ion of interest of rank 2
%          V4  - is a 9-component vector with the parameters of rank 4
%          V6  - is a 13-component vector with the parameters of rank 6
%                The format for each vector is from V_{-(rank)} to V_{+(rank)}
%
% Outputs: Hcf - A (2J+1)x(2J+1) matrix representing the Crystal Field Hamiltonian
%
% This function uses the normalised operators of Peter Fabi. Use the function
% norm_cfpars.m to convert Wybourne normalisation Bkq parameters to Vkq parameters
% used here. See FOCUS manual, Rutherford Appleton Laboratory report TR-95-023
%
% This function is based on the program ENGYLVL.BAS by Newman and Ng
% From their book: Crystal Field Handbook (CUP 2000).

% By Duc Le - Wed Aug 16 10:47:30 BST 2006 - duc.le@ucl.ac.uk
% mdl - updated 060821 - Non-zero component matrices were wrong... Oops! Now fixed.

% This file is part of the SAfiCF package. 
% Licenced under the GNU GPL v2 or later. 

% Checks that we have the right arguments.
if nargin == 2
  if ~iscell(V2)
    error(['If using only two arguments, V must be a cell array. Type "help ' ...
           mfilename '".']);
  else
    V4 = V2{2}; V6 = V2{3}; V2 = V2{1};
    if length(V2) ~= 5 || length(V4) ~= 9 || length(V6) ~=13
      error(['Invalid crystal field parameters. Type "help ' mfilename '".'])
    end
  end
elseif nargin == 4
  if ~isvector(V2) || length(V2) ~= 5
    error(['Invalid crystal field parameters V2. Type "help ' mfilename '".'])
  elseif ~isvector(V4) || length(V4) ~= 9
    error(['Invalid crystal field parameters V4. Type "help ' mfilename '".'])
  elseif ~isvector(V6) || length(V6) ~= 13
    error(['Invalid crystal field parameters V6. Type "help ' mfilename '".'])
  end
else
  error(['This function requires either four or two arguments. Type "help ' ...
         mfilename '".'])
end
if ~isscalar(J) || mod((J*2),1)
  error(['Total angular momentum J must be a scalar and integer or half integer.' ...
         ' Type "help ' mfilename '".'])
end

% The crystal field Hamiltonian operator in this normalisation is given by:
%
%                      k         q  k                   k                 k         q  k
%                     O    - (-1)  O                   O                 O    + (-1)  O      
%         ---   k      |q|          -|q|      ---  k    0       ---   k   |q|          -|q|  
% H   = i >    V     --------------------  +  >   V  ------  +  >    V  --------------------
%  cf     ---   -|q|         k                ---  0    k       ---   q           k
%        k,q<0            ||O   ||             k     ||O ||    k,q>0           ||O ||
%                            |q|                        0                         q
%
% where O_k^q are the Stevens operators of rank k.
%
% The normalisation is given by:
%                                                         2
%  k      k    k               k  2     ---    |<i|Okq|j>|
% V   =  B  ||O ||   where  ||O ||   =  >      ------------
%  q      q    q               q        ---i,j     2J+1
%
% Reference: FOCUS manual by Peter Fabi. Rutherford Appleton Laboratory Report RAL-TR-95-023 
% Reference: C. Rudowicz, J. Phys. C: Solid State Phys., vol 18, pp1415-1430 (1985).
%
% The energy matrix elements <LSJM_j|Okq|L'SJ'M_j'> are given by summing
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
% (j||O ||j) = 2   | -------------
%                 \|   (2J - k)!
%
% Reference: D.Smith and J.H.M. Thornley, Proc. Phys. Soc., 1966, vol 89, pp779.

% Implementation note: Because of the need to compute the normalisation constants, this 
% function would take a much longer time to complete than cf_hmltn.m so some ways of 
% reducing computation to speed up the function have been included. They maybe backported 
% into cf_hmltn.m in a future version - mdl 060816.

if 2*J-2>0
  RM2 = (1/4) * sqrt( factorial(2*J+2+1) / factorial(2*J-2) );
else % For J=0,0.5,1 CF Hamiltonian is zero.
  M = zeros(2*J+1);
  return;  
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
% Eg. Order 0 terms has it along the diagonal. Order 1 terms has it one element right,
% and so on. We can use this to calculate only for the orders with nonzero parameters.

% Zero order terms: B20, B40, B60
if sum(find(abs(find(V2)-3)==0)) | sum(find(abs(find(V4)-5)==0)) | sum(find(abs(find(V6)-7)==0))
  M20 = zeros(dimj); M40 = zeros(dimj); M60 = zeros(dimj);
  for i = 1:dimj
    Mj = i-1 - J;
    % Rank 2
    M20(i,i) = M20(i,i) + (-1)^(J-Mj) * threej([J 2 J; -Mj  0 Mj]) * RM2 * 2;
    % Rank 4
    M40(i,i) = M40(i,i) + (-1)^(J-Mj) * threej([J 4 J; -Mj  0 Mj]) * RM4 * 8; 
    % Rank 6
    M60(i,i) = M60(i,i) + (-1)^(J-Mj) * threej([J 6 J; -Mj  0 Mj]) * RM6 * 16;
  end
  norm2 = sqrt( trace(M20*conj(M20)) / (2*J+1) ); if norm2~=0; M20 = M20 .* V2(3) ./ norm2; end
  norm4 = sqrt( trace(M40*conj(M40)) / (2*J+1) ); if norm4~=0; M40 = M40 .* V4(5) ./ norm4; end
  norm6 = sqrt( trace(M60*conj(M60)) / (2*J+1) ); if norm6~=0; M60 = M60 .* V6(7) ./ norm6; end
  M = M + M20 + M40 + M60;
end % if nonzero zeroth order.

% One order terms: B21, B2-1, B41, B4-1, B61, B6-1
if sum(find(abs(find(V2)-3)==1)) | sum(find(abs(find(V4)-5)==1)) | sum(find(abs(find(V6)-7)==1))
  M21  = zeros(dimj); M41  = zeros(dimj); M61  = zeros(dimj);
  M2m1 = zeros(dimj); M4m1 = zeros(dimj); M6m1 = zeros(dimj);
  for i = 1:dimj-1 
    j = i + 1;
    Mj = i-1 - J;
    Mjp = i - J;
    % Calculates the threej values here to save computation time later.
    tj21 = threej([J 2 J; -Mj  1 Mjp]);     tj2m1 = threej([J 2 J; -Mj -1 Mjp]);
    tj41 = threej([J 4 J; -Mj  1 Mjp]);     tj4m1 = threej([J 4 J; -Mj -1 Mjp]);
    tj61 = threej([J 6 J; -Mj  1 Mjp]);     tj6m1 = threej([J 6 J; -Mj -1 Mjp]);
    % Rank 2
    M2m1(i,j) = M2m1(i,j) + (-1)^(J-Mj)   * tj21  * RM2 * -sqrt(-6);
    M2m1(i,j) = M2m1(i,j) - (-1)^(J-Mj-1) * tj2m1 * RM2 * -sqrt(-6);
    M21(i,j)  = M21(i,j)  + (-1)^(J-Mj)   * tj21  * RM2 * -sqrt(6);
    M21(i,j)  = M21(i,j)  + (-1)^(J-Mj+1) * tj2m1 * RM2 * -sqrt(6);
    % Rank 4
    M4m1(i,j) = M4m1(i,j) + (-1)^(J-Mj)   * tj41  * RM4 * -2/sqrt(-5);
    M4m1(i,j) = M4m1(i,j) - (-1)^(J-Mj-1) * tj4m1 * RM4 * -2/sqrt(-5);
    M41(i,j)  = M41(i,j)  + (-1)^(J-Mj)   * tj41  * RM4 * -2/sqrt(5);
    M41(i,j)  = M41(i,j)  + (-1)^(J-Mj+1) * tj4m1 * RM4 * -2/sqrt(5);
    % Rank 6
    M6m1(i,j) = M6m1(i,j) + (-1)^(J-Mj)   * tj61  * RM6 * -8/sqrt(-42);
    M6m1(i,j) = M6m1(i,j) - (-1)^(J-Mj-1) * tj6m1 * RM6 * -8/sqrt(-42);
    M61(i,j)  = M61(i,j)  + (-1)^(J-Mj)   * tj61  * RM6 * -8/sqrt(42);
    M61(i,j)  = M61(i,j)  + (-1)^(J-Mj+1) * tj6m1 * RM6 * -8/sqrt(42);
  end
  % Calculates the lower triangle of each matrix
  M21 = M21 + M21';    M41 = M41 + M41';    M61 = M61 + M61';
  M2m1 = M2m1 + M2m1'; M4m1 = M4m1 + M4m1'; M6m1 = M6m1 + M6m1';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M21*conj(M21)) / (2*J+1)); if norm~=0; M21 = M21 .* V2(4) ./ norm; end
  norm = sqrt(trace(M41*conj(M41)) / (2*J+1)); if norm~=0; M41 = M41 .* V4(6) ./ norm; end
  norm = sqrt(trace(M61*conj(M61)) / (2*J+1)); if norm~=0; M61 = M61 .* V6(8) ./ norm; end
  norm = sqrt(trace(M2m1*conj(M2m1))/(2*J+1)); if norm~=0; M2m1 = M2m1.*V2(2) ./ norm; end
  norm = sqrt(trace(M4m1*conj(M4m1))/(2*J+1)); if norm~=0; M4m1 = M4m1.*V4(4) ./ norm; end
  norm = sqrt(trace(M6m1*conj(M6m1))/(2*J+1)); if norm~=0; M6m1 = M6m1.*V6(6) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M21 + M41 + M61 + sqrt(-1).*(M2m1 + M4m1 + M6m1);
end % if nonzero one order.

% Two order terms: B22, B2-2, B42, B4-2, B62, B6-2
if sum(find(abs(find(V2)-3)==2)) | sum(find(abs(find(V4)-5)==2)) | sum(find(abs(find(V6)-7)==2))
  M22  = zeros(dimj); M42  = zeros(dimj); M62  = zeros(dimj);
  M2m2 = zeros(dimj); M4m2 = zeros(dimj); M6m2 = zeros(dimj);
  for i = 1:dimj-2 
    j = i + 2;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj22 = threej([J 2 J; -Mj  2 Mjp]);     tj2m2 = threej([J 2 J; -Mj -2 Mjp]);
    tj42 = threej([J 4 J; -Mj  2 Mjp]);     tj4m2 = threej([J 4 J; -Mj -2 Mjp]);
    tj62 = threej([J 6 J; -Mj  2 Mjp]);     tj6m2 = threej([J 6 J; -Mj -2 Mjp]);
    % Rank 2
    M2m2(i,j) = M2m2(i,j) + (-1)^(J-Mj)   * tj22  * RM2 * 2/sqrt(-6);
    M2m2(i,j) = M2m2(i,j) - (-1)^(J-Mj-2) * tj2m2 * RM2 * 2/sqrt(-6);
    M22(i,j)  = M22(i,j)  + (-1)^(J-Mj)   * tj22  * RM2 * 2/sqrt(6);
    M22(i,j)  = M22(i,j)  + (-1)^(J-Mj+2) * tj2m2 * RM2 * 2/sqrt(6);
    % Rank 4
    M4m2(i,j) = M4m2(i,j) + (-1)^(J-Mj)   * tj42  * RM4 * 4/sqrt(-10);
    M4m2(i,j) = M4m2(i,j) - (-1)^(J-Mj-2) * tj4m2 * RM4 * 4/sqrt(-10);
    M42(i,j)  = M42(i,j)  + (-1)^(J-Mj)   * tj42  * RM4 * 4/sqrt(10);
    M42(i,j)  = M42(i,j)  + (-1)^(J-Mj+2) * tj4m2 * RM4 * 4/sqrt(10);
    % Rank 6
    M6m2(i,j) = M6m2(i,j) + (-1)^(J-Mj)   * tj62  * RM6 * 16/sqrt(-105);
    M6m2(i,j) = M6m2(i,j) - (-1)^(J-Mj-2) * tj6m2 * RM6 * 16/sqrt(-105);
    M62(i,j)  = M62(i,j)  + (-1)^(J-Mj)   * tj62  * RM6 * 16/sqrt(105);
    M62(i,j)  = M62(i,j)  + (-1)^(J-Mj+2) * tj6m2 * RM6 * 16/sqrt(105);
  end
  % Calculates the lower triangle of each matrix
  M22 = M22 + M22';    M42 = M42 + M42';    M62 = M62 + M62';
  M2m2 = M2m2 + M2m2'; M4m2 = M4m2 + M4m2'; M6m2 = M6m2 + M6m2';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M22*conj(M22)) / (2*J+1)); if norm~=0; M22 = M22 .* V2(5) ./ norm; end
  norm = sqrt(trace(M42*conj(M42)) / (2*J+1)); if norm~=0; M42 = M42 .* V4(7) ./ norm; end
  norm = sqrt(trace(M62*conj(M62)) / (2*J+1)); if norm~=0; M62 = M62 .* V6(9) ./ norm; end
  norm = sqrt(trace(M2m2*conj(M2m2))/(2*J+1)); if norm~=0; M2m2 = M2m2.*V2(1) ./ norm; end
  norm = sqrt(trace(M4m2*conj(M4m2))/(2*J+1)); if norm~=0; M4m2 = M4m2.*V4(3) ./ norm; end
  norm = sqrt(trace(M6m2*conj(M6m2))/(2*J+1)); if norm~=0; M6m2 = M6m2.*V6(5) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M22 + M42 + M62 + sqrt(-1).*(M2m2 + M4m2 + M6m2);
end % if nonzero two order

% Three order terms: B43, B4-3, B63, B6-3
if sum(find(abs(find(V4)-5)==3)) | sum(find(abs(find(V6)-7)==3))
  M43  = zeros(dimj); M63  = zeros(dimj);
  M4m3 = zeros(dimj); M6m3 = zeros(dimj);
  for i = 1:dimj-3 
    j = i + 3;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj43 = threej([J 4 J; -Mj  3 Mjp]);     tj4m3 = threej([J 4 J; -Mj -3 Mjp]);
    tj63 = threej([J 6 J; -Mj  3 Mjp]);     tj6m3 = threej([J 6 J; -Mj -3 Mjp]);
    % Rank 4
    M4m3(i,j) = M4m3(i,j) + (-1)^(J-Mj)   * tj43  * RM4 * -2/sqrt(-35);
    M4m3(i,j) = M4m3(i,j) - (-1)^(J-Mj-3) * tj4m3 * RM4 * -2/sqrt(-35);
    M43(i,j)  = M43(i,j)  + (-1)^(J-Mj)   * tj43  * RM4 * -2/sqrt(35);
    M43(i,j)  = M43(i,j)  + (-1)^(J-Mj+3) * tj4m3 * RM4 * -2/sqrt(35);
    % Rank 6
    M6m3(i,j) = M6m3(i,j) + (-1)^(J-Mj)   * tj63  * RM6 * -8/sqrt(-105);
    M6m3(i,j) = M6m3(i,j) - (-1)^(J-Mj-3) * tj6m3 * RM6 * -8/sqrt(-105);
    M63(i,j)  = M63(i,j)  + (-1)^(J-Mj)   * tj63  * RM6 * -8/sqrt(105);
    M63(i,j)  = M63(i,j)  + (-1)^(J-Mj+3) * tj6m3 * RM6 * -8/sqrt(105);
  end
  % Calculates the lower triangle of each matrix
  M43 = M43 + M43';    M63 = M63 + M63';
  M4m3 = M4m3 + M4m3'; M6m3 = M6m3 + M6m3';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M43*conj(M43)) / (2*J+1)); if norm~=0; M43 = M43 .* V4(8) ./ norm; end
  norm = sqrt(trace(M63*conj(M63)) / (2*J+1)); if norm~=0; M63 = M63 .* V6(10) ./ norm; end
  norm = sqrt(trace(M4m3*conj(M4m3))/(2*J+1)); if norm~=0; M4m3 = M4m3.*V4(2) ./ norm; end
  norm = sqrt(trace(M6m3*conj(M6m3))/(2*J+1)); if norm~=0; M6m3 = M6m3.*V6(4) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M43 + M63 + sqrt(-1).*(M4m3 + M6m3);
end % if nonzero three order

% Four order terms: B44, B4-4, B64, B6-4
if sum(find(abs(find(V4)-5)==4)) | sum(find(abs(find(V6)-7)==4))
  M44  = zeros(dimj); M64  = zeros(dimj);
  M4m4 = zeros(dimj); M6m4 = zeros(dimj);
  for i = 1:dimj-4 
    j = i + 4;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj44 = threej([J 4 J; -Mj  4 Mjp]);     tj4m4 = threej([J 4 J; -Mj -4 Mjp]);
    tj64 = threej([J 6 J; -Mj  4 Mjp]);     tj6m4 = threej([J 6 J; -Mj -4 Mjp]);
    % Rank 4
    M4m4(i,j) = M4m4(i,j) + (-1)^(J-Mj)   * tj44  * RM4 * 8/sqrt(-70);
    M4m4(i,j) = M4m4(i,j) - (-1)^(J-Mj-4) * tj4m4 * RM4 * 8/sqrt(-70);
    M44(i,j)  = M44(i,j)  + (-1)^(J-Mj)   * tj44  * RM4 * 8/sqrt(70);
    M44(i,j)  = M44(i,j)  + (-1)^(J-Mj+4) * tj4m4 * RM4 * 8/sqrt(70);
    % Rank 6
    M6m4(i,j) = M6m4(i,j) + (-1)^(J-Mj)   * tj64  * RM6 * 16/3/sqrt(-14);
    M6m4(i,j) = M6m4(i,j) - (-1)^(J-Mj-4) * tj6m4 * RM6 * 16/3/sqrt(-14);
    M64(i,j)  = M64(i,j)  + (-1)^(J-Mj)   * tj64  * RM6 * 16/3/sqrt(14);
    M64(i,j)  = M64(i,j)  + (-1)^(J-Mj+4) * tj6m4 * RM6 * 16/3/sqrt(14);
  end
  % Calculates the lower triangle of each matrix
  M44 = M44 + M44';    M64 = M64 + M64';
  M4m4 = M4m4 + M4m4'; M6m4 = M6m4 + M6m4';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M44*conj(M44)) / (2*J+1)); if norm~=0; M44 = M44 .* V4(9) ./ norm; end
  norm = sqrt(trace(M64*conj(M64)) / (2*J+1)); if norm~=0; M64 = M64 .* V6(11)./ norm; end
  norm = sqrt(trace(M4m4*conj(M4m4))/(2*J+1)); if norm~=0; M4m4 = M4m4.*V4(1) ./ norm; end
  norm = sqrt(trace(M6m4*conj(M6m4))/(2*J+1)); if norm~=0; M6m4 = M6m4.*V6(3) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M44 + M64 + sqrt(-1).*(M4m4 + M6m4);
end % if nonzero four order.

% Five order terms: B65, B6-5
if sum(find(abs(find(V6)-7)==5))
  M65  = zeros(dimj);
  M6m5 = zeros(dimj);
  for i = 1:dimj-5 
    j = i + 5;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj65 = threej([J 6 J; -Mj  5 Mjp]);     tj6m5 = threej([J 6 J; -Mj -5 Mjp]);
    % Rank 6
    M6m5(i,j) = M6m5(i,j) + (-1)^(J-Mj)   * tj65  * RM6 * -8/3/sqrt(-77);
    M6m5(i,j) = M6m5(i,j) - (-1)^(J-Mj-5) * tj6m5 * RM6 * -8/3/sqrt(-77);
    M65(i,j)  = M65(i,j)  + (-1)^(J-Mj)   * tj65  * RM6 * -8/3/sqrt(77);
    M65(i,j)  = M65(i,j)  + (-1)^(J-Mj+5) * tj6m5 * RM6 * -8/3/sqrt(77);
  end
  % Calculates the lower triangle of each matrix
  M65 = M65 + M65'; M6m5 = M6m5 + M6m5';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M65*conj(M65)) / (2*J+1)); if norm~=0; M65 = M65 .* V6(12)./ norm; end
  norm = sqrt(trace(M6m5*conj(M6m5))/(2*J+1)); if norm~=0; M6m5 = M6m5.*V6(2) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M65 + sqrt(-1).*M6m5;
end % if nonzero five order.

% Six order terms: B66, B6-6
if sum(find(abs(find(V6)-7)==6))
  M66  = zeros(dimj);
  M6m6 = zeros(dimj);
  for i = 1:dimj-6
    j = i + 6;
    Mj = i-1 - J;
    Mjp = j-1 - J;
    % Calculates the threej values here to save computation time later.
    tj66 = threej([J 6 J; -Mj  6 Mjp]);     tj6m6 = threej([J 6 J; -Mj -6 Mjp]);
    % Rank 6
    M6m6(i,j) = M6m6(i,j) + (-1)^(J-Mj)   * tj66  * RM6 * 16/sqrt(-231);
    M6m6(i,j) = M6m6(i,j) - (-1)^(J-Mj-6) * tj6m6 * RM6 * 16/sqrt(-231);
    M66(i,j)  = M66(i,j)  + (-1)^(J-Mj)   * tj66  * RM6 * 16/sqrt(231);
    M66(i,j)  = M66(i,j)  + (-1)^(J-Mj+6) * tj6m6 * RM6 * 16/sqrt(231);
  end
  % Calculates the lower triangle of each matrix
  M66 = M66 + M66'; M6m6 = M6m6 + M6m6';
  % Calculates the normalisation constant and apply it to the matrix
  norm = sqrt(trace(M66*conj(M66)) / (2*J+1)); if norm~=0; M66 = M66 .* V6(13) ./ norm; end
  norm = sqrt(trace(M6m6*conj(M6m6)) / (2*J+1)); if norm~=0; M6m6 = M6m6 .* V6(1) ./ norm; end
  % The normalisation for negative q is imaginary so we need to add i back in explicitly
  M = M + M66 + sqrt(-1).*M6m6;
end % if nonzero six order
