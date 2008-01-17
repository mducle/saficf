function out = fvoigt(xdat,pars)
% Calculates a pseudo-voigt function with given parameters at given ordinate.
%
% Syntax:  out  = fvoigt(xdat,pars)
%
% Inputs:  xdat - vector - is the ordinate at which to calculate
%          pars - vector - [centre fwhm height lfrac] - lfrac is lorentzian fraction
%
% Outputs: out  - vector - is the abscisa which is calculated.

% Reference: pvoigt.m file from the fitting routines for ID20 by SPC, SBW.

% Makes equation looks better
c = pars(1);
w = pars(2);
h = pars(3);
f = pars(4);
x = xdat(:);

out  = h .* ( f./(1+4*((x-c)/w).^2) + (1-f).*exp(-4.*log(2).*((x-c)/w).^2) );
