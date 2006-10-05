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

if ~exist('maxsplit')  maxsplit = 50;   end
if ~exist('J')         J        = 4;    end

% Gets the allowed parameters
allowed = ptgp(ptgpstr);

% Checks that the site symmetry is cubic. If so, need  V44 = 5.V40; V64 = 21.V60;
if strcmp(ptgpstr,'cubic') | strcmp(lower(ptgpstr),'t') | strcmp(lower(ptgpstr),'o')
  B = rand(1,2) .* splitfactor;
  Vstart = {zeros(1,5) [0 0 0 0 B(1) 0 0 0 5*B(1)] [zeros(1,6) B(2) 0 0 0 21*B(2) 0 0]}; 
else
  for ind_sites = 1:size(allowed,2)
    % Works out the maximum splitting
    E = eig(norm_cfhmltn(J,allowed{:,ind_sites})); 
    split_unity = max(E'-min(E));         % Full CF split with all parameters set to unity.
    splitfactor = maxsplit / split_unity; % If all pars = splitfactor, full split = maxsplit
    for ind_k = 1:3
      % The splitfactor ensures that the energy split is less than maxsplit
      Vstart{ind_k,ind_sites} = rand(1,4*ind_k+1) .* splitfactor .* allowed{ind_k,ind_sites};
    end
  end
end
