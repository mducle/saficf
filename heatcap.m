function [C, J] = heatcap(En,T);
% heatcap - Calculates the heat capacity as function of temperature
%
% Syntax:  [C, J] = heatcap(En,T)
% 
% Inputs:  En is a vector with the energy levels in meV;
%          T is a vector containing the temperature at which to sample.
% Outputs: C is the heat capacity as a function of temperature C(T);
%          J is the Jacobian: Jij = dC(Ti)/dEj 

% Uses the thermodynamic relation:                            2
%     dU|           1   d    1 dZ      1       [ 1 ---    -BEn ]   1 ---  2 -BEn 
% C = --|        = ---- -- [ - -- ] = ---- { - [ - >   E e     ] + - >   E e     }
%     dT|N,{V,p}   kT^2 dB   Z dB     kT^2     [ Z ---n n      ]   Z ---n n
%                   ^          ^                    ^                 ^
% where B = 1/kT    ^---dB/dT  ^---U        <<U>> --^       <<U^2>> --^
% where Z is the partition function [sum(exp(-En/kT))] and U is the internal
% energy given by [sum(En * exp(-En/kT)) / Z].

% By Duc Le 2005 - duc.le@ucl.ac.uk
% mdl - updated 070129 - minor fixes: size->length(T)

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
k_B = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e = 1.60217653e-19;    % C - Charge of electron
N_A = 6.0221415e23;      % Avogadro's number

% Sorts out dimensionality problems
T = T(:)';

% Converts energy levels from meV to J
En = En-min(En); En = En .* (Q_e/1000);

% Some definitions
beta = 1 ./ (k_B*T);

for ind = 1:length(T)
% Finds the partition function Z(T)
  Z(ind) = sum( exp(-beta(ind) .* En) );  % dimensionless

% Finds the expectation of the internal energy: <<U(T)>> 
  U(ind) = sum( En .* exp(-beta(ind) .* En) ) ./ Z(ind);  % in J

% Finds the expectation of the square of the  energy <<U^2(T)>>
  U2(ind) = sum( En.^2 .* exp(-beta(ind) .* En) ) ./ Z(ind);  % in J^2
end

C = (U2 - U.^2) ./ (k_B .* T.^2);  % in J/K/atom

C = C' .* N_A;  % Convert to J/K/mol

% Checks if we have to calculate the Jacobian.
if nargout > 1
% The Jacobian is given by:
%
%      dC(Ti)    1   1   [    -BiEj    2   -BiEj ]      [ -2BiEj      -2BiEj ]
% J  = ------ = ---- - { [2E e      - E B e      ] - 2E [e       - E e       ] }
%  ij   dEj     kT^2 Z   [  j          j i       ]     j[           j        ]
%
% The terms in square brackets are dU/dE and dU2/dE respectively, and will be
% denoted by dU and dU2 in the equations below for clarity.

% Converts En vector into row form
  if size(En,1) ~=1
    En = En';
  end
  for ind = 1:size(T,2)
    dU(ind,:) = 2.*En.*exp(-beta(ind).*En) - En.^2.*beta(ind).*exp(-beta(ind).*En);
    dU2(ind,:) = exp(-2.*beta(ind).*En) - En.*exp(-2.*beta(ind).*En);
    J(ind,:) = ( beta(ind)/T(ind) / Z(ind) ) .* ( dU(ind,:) - 2.*En.*dU2(ind,:) );
  end
else  % Don't work out the Jacobian
  J = []; 
end





