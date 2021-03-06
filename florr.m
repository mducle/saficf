function out = florr(xdat,pars)
% Calculates a lorentzian function with given parameters at given ordinate.
%
% Syntax:  out = florr(pars,xdat)
% Inputs:  xdat - vector - is the ordinate at which to calculate
%          pars - vector - [centre fwhm area]
%
% Outputs: out  - vector - is the abscisa which is calculated.

% Makes equation looks better
c = pars(1);
w = pars(2);
a = pars(3);
x = xdat(:);

%out = (h*(w/2)^2) ./ ( (x-c).^2+(w/2)^2 );
%out = (2*a/w/pi) ./ (1 + 4*((x-c)/w).^2);

out = (a/pi*(w/2)) ./ ( (x-c).^2+(w/2)^2 );
