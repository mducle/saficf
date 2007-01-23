function range = saficf_range(ptgpstr,maxsplit,J)
% Determines the range on the normalised CF parameters given a maximum splitting.
%
% Syntax:  range = saficf_range(ptgpstr,maxsplit,J)
%
% Inputs:  ptgpstr  - string - point group symmetry of the magnetic ion of interest
%          maxsplit - scalar - energy of the highest CF state (meV), optional
%          J        - scalar - total angular momentum quantum number, optional
%
% Outputs: range    - vector - num_sites sized vector of ranges [-range,+range]
%
% If not specified, maxsplits defaults to 50meV, and J to 4.

% mdl - update 070123 - Fixed cubic Fabi normalised parameters.

if ~exist('maxsplit')  maxsplit = 50;   end
if ~exist('J')         J        = 4;    end

% Gets the allowed parameters
allowed = ptgp(ptgpstr);

for ind_sites = 1:size(allowed,2)
  % Works out the maximum splitting
  % Checks that the site symmetry is cubic. If so, need  V44 = 5.V40; V64 = -21.V60;
  if strcmp(ptgpstr,'cubic') | strncmp(lower(ptgpstr),'t',1) | strncmp(lower(ptgpstr),'o',1)
    if (J<2.5)        % Means E = zeros(1,2*J+1), so split_unity = 0, get Inf error.
      maxsplit = 0;
      E = 1:(2*J+1);
    else
      V = norm_cfpars(J,{zeros(1,5) [0 0 0 0 1 0 0 0 5] [zeros(1,6) 1 0 0 0 -21 0 0]}); 
      E = eig(norm_cfhmltn(J,V));
    end
  else
    E = eig(norm_cfhmltn(J,allowed{:,ind_sites})); 
  end
  split_unity = max(E'-min(E));              % Full CF split with all parameters set to unity.
  range(ind_sites) = maxsplit / split_unity; % If all pars = splitfactor, full split = maxsplit
end
