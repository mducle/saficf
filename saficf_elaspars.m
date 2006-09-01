function pars = saficf_elaspars(xdat,ydat)
% Estimates the centre, FWHM and height of the elastic peaks, and best fitted lineshape
%
% Syntax:  pars = saficf_elaspars(xdat,ydat)
%
% Inputs:  xdat - vector     - the input independent variable (energy transfer in meV)
%          ydat - vector     - the input dependent variable (intensity or S(q,w))
%
% Outputs: pars - cell array - {type_str [centre fwhm height lfrac]} - parameters of fit.
%                              type_str is a string: 'gauss', 'lor' or 'pvoigt'
%                              pars{2} is a vector of parameters suitable to lineshape.
%
% This function uses the Levenberg-Marquardt algorithm implemented by Shrager, Jutan and
%   Muzic in the file speclsqr.m from the spec1d package to test for the lineshape that 
%   best matches the elastic peak, from choice of Gaussian, Lorentzian, and a pseudo-Voigt
%   function (superposition of Gaussian and Lorentzian) with the centres and height fixed,
%   with lfrac as the ratio of Lorentzian to Gaussian FWHM.

% Determine the height and width of the elastic peak. 
elas_range = ydat(find(abs(xdat)<5));            % Only consider Et from -5:5 meV.
elas_height = max(elas_range);                   % Finds max value of S(q,w) in range.
elas_centre = xdat(find(ydat==elas_height));     % Finds the Et corresponding to peak.
elas_tophalf = elas_range(find(elas_range>(elas_height/2.5))); % All S(q,w) > half max
elas_tophalf = sort(elas_tophalf);
% Finds the FWHM by adding the Et of the values of S(q,w) near half maxima.
elas_fwhm = sum(abs([xdat(ydat==elas_tophalf(1)) xdat(ydat==elas_tophalf(2))]));
% Puts it all into a vector of parameters for a pseudo-voigt and gaussian function
voigt_p = [elas_centre elas_fwhm elas_height*elas_fwhm 0.1];
gauss_p = voigt_p(1:3);
% So that speclsqr doesn't complain
xdat = xdat(:); ydat = ydat(:);

% Only fits up to FWHM on the energy loss side to ignore quasi-elastic scattering.
ind_whi = max([find(ydat==elas_tophalf(1)) find(ydat==elas_tophalf(2))]);
x = xdat(1:ind_whi);
y = ydat(1:ind_whi);

% Fits the parameters using Levenberg-Marquardt least squares fitting.
pars_voigt = speclsqr(x,y,ones(1,length(y)),voigt_p,ones(1,4)*0.1,@fvoigt);
pars_gauss = speclsqr(x,y,ones(1,length(y)),gauss_p,ones(1,3)*0.1,@fgauss);
pars_lorr  = speclsqr(x,y,ones(1,length(y)),gauss_p,ones(1,3)*0.1,@florr );

% Determine the best fitted lineshape by calculating the least squares difference
lsq(1) = sum( (feval(@fvoigt,x,pars_voigt) - y).^2 );
lsq(2) = sum( (feval(@fgauss,x,pars_gauss) - y).^2 );
lsq(3) = sum( (feval(@florr ,x,pars_gauss) - y).^2 );
switch find(lsq == min(lsq))
  case 1; pars = {'fvoigt',pars_voigt};
  case 2; pars = {'fgauss',pars_gauss};
  case 3; pars = {'florr' ,pars_gauss};
end

