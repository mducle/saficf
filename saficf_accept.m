function aflag = saficf_accept(cost,c)
% Tests whether to accept a new configuration or not based on FSA criteria
%
% Syntax:  aflag = saficf_accept(cost,c)
%
% Inputs:  cost  - scalar - the cost difference between the new and old configurations
%          c     - scalar - the control parameter (temperature)
%
% Outputs: aflag - logical - 1=accept, 0=reject

if cost < 0                          % Accepts all transitions that lowers the cost
  aflag = 1;
elseif rand < (1/(1 + exp(cost/c)))  % Accepts transitions with a Fermi-Dirac probability
  aflag = 1;
else                                 % Rejects all other transitions
  aflag = 0;
end
