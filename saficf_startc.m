function cstart = saficf_startc(xdat,ydat,edat,J,T,Ei,freq,ptgpstr,elas_peak)
% Defines the start temperature so that at least 49% of transitions are accepted.
%
% Syntax:  cstart = saficf_startc(xdat,ydat,J,T,Ei,freq,ptgpstr,elas_peak)
% 
% Inputs:  xdat      - vector - input independent variable (energy transfer in meV)
%          ydat      - vector - input dependent variable (intensity or S(q,w))
%          edat      - vector - input error on y-data. 
%          J         - scalar - total angular momentum quantum number
%          T         - scalar - temperature of sample (Kelvin)
%          Ei        - scalar - incident energy (meV)
%          freq      - scalar - chopper frequency (Hz)
%          ptgpstr   - string - point group of magnetic ion of interest.
%          elas_peak - vector - the elastic peak evaluated at values of xdat
%
% Outputs: cstart  - scalar - the start control parameter value

% The acceptance criteria used by the fast annealing algorithm is:
%                1
% P     =  -------------   Ref: Szu and Harley, Phys. Lett. A 122 (1987) 157, eqn 15
%  FSA     1 + exp(dE/c)
%
% So to ensure 49% acceptance, we need c_start = dE / log((1-n)/n), where n = 0.49

% Generates and finds the cost of some random configurations.
%   NB. Each loop takes approximately 0.12s on a 1.6Ghz Pentium M for triclinic
%       and approximately 0.05s for hexagonal symmetry.

num_set = length(T);

for i_set = 1:num_set
  if max(max(abs(edat{i_set}))) < 1e-5
    costflag = 1;                                 % Cost is root mean square difference
  else
    costflag = 0;                                 % Cost is chi-square
  end
end

for ind_conf = 1:10

  cost(ind_conf) = 0;
  V = saficf_genstart(ptgpstr,50,J);              % Generates new configuration
  for i_set = 1:num_set                           %   and new intensity factor
    intfac{i_set} = rand * max(ydat{i_set}) / 10;
    spectmp = saficf_genspec(J,T(i_set),V,xdat{i_set},Ei(i_set),freq(i_set),'fgauss');
    spec = elas_peak{i_set} ...                   % Add elastic and inelastic peaks
           + spectmp.*intfac{i_set};
    if costflag                                   % Calculates new cost 
      cost(ind_conf) = cost(ind_conf) + sqrt( sum( (spec - ydat{i_set}).^2 ) );
    else
      cost(ind_conf) = cost(ind_conf) + sum( (spec - ydat{i_set}).^2 ./ (edat{i_set}.^2) );
    end
  end


%  V = saficf_genstart(ptgpstr,50,J);
%  intfac = rand*max(ydat)/10;
%  spec = saficf_genspec(J,T,V,xdat(:)',Ei,freq,'fgauss') .* intfac + elas_peak(:)';
%  cost(ind_conf) = sqrt( sum( (spec - ydat(:)').^2 ) );

%  specs{ind_conf} = spec;
end

mean_cost = sum(cost)/10;

cstart = mean_cost / 25;   % 1/log((1-n)/n) = 25 for n = 0.49
%cstart = specs;
