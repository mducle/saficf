function Vnew = saficf_perturb(V,c,range)
% Generates new state from the set of CF parameters V.
%
% Syntax:  Vnew  = saficf_perturb(V,c)
%
% Inputs:  V     - cell array - {V2 V4 V6} of Fabi normalised CF parameters
%          c     - scalar     - the control parameter (temperature).
%          range - vector     - the range of CF parameters on each site.
%
% Outputs: Vnew  - cell array - {Vn2 Vn4 Vn6} of the new Cauchy distributed parameters

% This function uses the criteria of the Fast Simulated Annealing algorithm, so 
%   the new parameters are Cauchy distributed.
% Reference: H. Szu and R. Hartley, Physics Letters A, 122 (1987) 157

% By Duc Le - Tue Aug 22 22:15:37 BST 2006 - duc.le@ucl.ac.uk

% mdl - updated 070123 - Fixed handling of cubic Fabi normalised parameters

for ind_sites = 1:size(V,2)
  Vnew(:,ind_sites) = {zeros(1,5); zeros(1,9); zeros(1,13)};

  % Checks that the site symmetry is cubic. If so, need  V44 = sqrt(5/7).V40; V64 = -sqrt(7).V60;
  %                                              That is B44 = 5.B40;         B64 = -21.B60;
  if   abs(V{2,ind_sites}(9)/V{2,ind_sites}(5)-sqrt(5/7)) < 1e-10 ...
    && abs(V{3,ind_sites}(11)/V{3,ind_sites}(7)+sqrt(7))  < 1e-10
    B = [V{2,ind_sites}(5)+lrnd(c) V{3,ind_sites}(7)+lrnd(c)];
    Vnew(:,ind_sites) = {zeros(1,5); [0 0 0 0 B(1) 0 0 0 sqrt(5/7)*B(1)]; ...
                         [zeros(1,6) B(2) 0 0 0 -sqrt(7)*B(2) 0 0]}; 
  else
    for ind_k = 1:3                             % Loop over the orders 2, 4, 6
      for ind_q = find(V{ind_k,ind_sites}(:))'  % Loop over nonzero components Bkq
        Vtmp = V{ind_k,ind_sites}(ind_q) + lrnd(c);
        while Vtmp > range(ind_sites) 
          Vtmp = V{ind_k,ind_sites}(ind_q) + lrnd(c);
        end
        Vnew{ind_k,ind_sites}(ind_q) = Vtmp; 
      end
    end
  end

end % ind_sites
