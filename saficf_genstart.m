function Vstart = saficf_genstart(ptgpstr,maxsplit,J)
% Generates a random set of crystal field parameters of appropriate symmetry
%
% Syntax:  Vstart = saficf_genstart(ptgpstr,maxsplit,J)
%
% Inputs:  ptgpstr  - string     - point group symmetry of the magnetic ion of interest
%          maxsplit - scalar     - energy of the highest CF state (meV), optional
%          J        - scalar     - total angular momentum quantum number, optional
%
% Outputs: Vstart   - cell array - {V2 V4 V6} randomly generated CF parameters
%
% If not specified, maxsplits defaults to 50meV, and J to 4.

% mdl - updated 070123 - fixed cubic parameters under Fabi normalisation.
%                      - fixed Vstart so negative CF parameters now generated.

if ~exist('maxsplit')  maxsplit = 50;   end
if ~exist('J')         J        = 4;    end

% Gets the allowed parameters
allowed = ptgp(ptgpstr);

for ind_sites = 1:size(allowed,2)
  % Checks that the site symmetry is cubic. If so, need  B44 = 5.B40; B64 = -21.B60;
  if strcmp(ptgpstr,'cubic') || strncmp(lower(ptgpstr),'t',1) || strncmp(lower(ptgpstr),'o',1)
    E = eig(cf_hmltn(J,zeros(1,5),[0 0 0 0 1 0 0 0 5],[zeros(1,6) 1 0 0 0 -21 0 0]));
    split_unity = max(E-min(E));            % Full CF split with all parameters set to unity.
    if (J<2.5)
      splitfactor = 0;                      % J<2.5: split_unity=0, so this avoids a Inf error.
    else
      splitfactor = maxsplit / split_unity; % If all pars = splitfactor, full split = maxsplit
    end
    B = rand(1,2) .* splitfactor;
    Vstart(:,ind_sites) = norm_cfpars(J,{zeros(1,5) [0 0 0 0 B(1) 0 0 0 5*B(1)] ...
                                         [zeros(1,6) B(2) 0 0 0 -21*B(2) 0 0]})'; 
  else
    % Works out the maximum splitting
    E = eig(norm_cfhmltn(J,allowed{:,ind_sites})); 
    split_unity = max(E-min(E));            % Full CF split with all parameters set to unity.
    splitfactor = maxsplit / split_unity;   % If all pars = splitfactor, full split = maxsplit
    for ind_k = 1:3
      % The splitfactor ensures that the energy split is less than maxsplit
      Vstart{ind_k,ind_sites} = (rand(1,4*ind_k+1).*2-1) .* splitfactor .* allowed{ind_k,ind_sites};
    end
  end
end
