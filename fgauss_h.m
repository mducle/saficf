function out = fgauss(xdat,pars)
% Calculates a gaussian function with given parameters at given ordinate.
%
% Syntax:  out = fgauss(pars,xdat)
%
% Inputs:  xdat - vector - is the ordinate at which to calculate
%          pars - vector - [centre fwhm height]
%
% Outputs: out  - vector - is the abscisa which is calculated.

% Makes equation looks better
c = pars(1);
w = pars(2);
h = pars(3);
x = xdat(:);

out = h .* exp(-4.*log(2).*((x-c)./w).^2);

