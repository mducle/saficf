function Tnew = saficf_schedul(Tstart,iteration)
% Generates a new value for the control parameter using the fast annealing schedule
%
% Syntax:  Tnew = saficf_schedul(Tstart,iteration)
%
% Inputs:  Tstart    - scalar - initial temperature
%          iteration - scalar - number of iterations so far
%
% Outputs: Tnew      - scalar - the new temperature

% Reference: H. Szu and R. Hartley, Physics Letters A, 122 (1987) 157

% By Duc Le 2006 - duc.le@ucl.ac.uk

Tnew = Tstart / (1 + iteration);
