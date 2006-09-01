function Vnew = saficf_perturb(V,c)
% Generates new state from the set of CF parameters V.
%
% Syntax:  Vnew = saficf_perturb(V,c)
%
% Inputs:  V    - cell array - {V2 V4 V6} of Fabi normalised CF parameters
%          c    - scalar     - the control parameter (temperature).
%
% Outputs: Vnew - cell array - {Vn2 Vn4 Vn6} of the new Cauchy distributed parameters

% This function uses the criteria of the Fast Simulated Annealing algorithm, so 
%   the new parameters are Cauchy distributed.
% Reference: H. Szu and R. Hartley, Physics Letters A, 122 (1987) 157

% By Duc Le - Tue Aug 22 22:15:37 BST 2006 - duc.le@ucl.ac.uk

Vnew = {zeros(1,5) zeros(1,9) zeros(1,13)};

% Checks that the site symmetry is cubic. If so, need  V44 = 5.V40; V64 = 21.V60;
if abs(V{2}(9)/V{2}(5)-5) < 1e-10 && abs(V{3}(11)/V{3}(7)-21) < 1e-10
  B = [V{2}(5)+lrnd(c) V{3}(7)+trnd(c)];
  Vnew = {zeros(1,5) [0 0 0 0 B(1) 0 0 0 5*B(1)] [zeros(1,6) B(2) 0 0 0 21*B(2) 0 0]}; 
else
  for ind_k = 1:3                % Loop over the orders 2, 4, 6
    for ind_q = find(V{ind_k})   % Loop over nonzero components Bkq
      Vnew{ind_k}(ind_q) = V{ind_k}(ind_q) + lrnd(c); 
    end
  end
end

