function stopflag = saficf_stopc(dEold,dEnew,stopflag)
% Determines when the SA temperature sequence should stop.
%   Criteria used is when energy has not changed for 5 previous control parameter steps 
%
% Syntax:  stopflag = saficf_stopc(dEold,dEnw,stopflag)
%
% Inputs:  dEold    - scalar - the old difference in energy
%          dEnew    - scalar - the new difference in energy
%          stopflag - scalar - the previous value of stop flag
%
% Outputs: stopflag - scalar - the new value of stop flag

if abs(dEold-dEnew) < 1e-10
  stopflag = stopflag + 1;
end
