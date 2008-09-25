function out = florr_h(xdat,pars)
% Calculates a lorentzian function with given parameters at given ordinate.
%
% Syntax:  out = florr(pars,xdat)
% Inputs:  xdat - vector - is the ordinate at which to calculate
%          pars - vector - [centre fwhm height]
%
% Outputs: out  - vector - is the abscisa which is calculated.

% Makes equation looks better
c = pars(1);
w = pars(2);
h = pars(3);
x = xdat(:);

out = (h*(w/2)^2) ./ ( (x-c).^2+(w/2)^2 );

