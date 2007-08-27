function [spec,fitpars,total_area,elas_area] = saficf_arbIrenorm(VHcf,xdat,ydat,edat,T,Ei,total_area,elas_area)
% Fits a set of CF parameters / a CF Hamiltonian to a spectra with arbitrary unit intensity
%
% Syntax:  [spec,pars,Atot,Ael] = saficf_arbIrenorm(VHcf,xdat,ydat,edat,Atot,Ael)
%
% Inputs:  VHcf - Either - square matrix = the CF hamiltonian
%                        - cell array = {J, V2, V4, V6} - the total angular momenturm q. num 
%                            and CF parameters in Stevens' normalisation
%                        - 2-column vector = a list of peaks / transition matrix elements pairs.
%          xdat - vector - The spectra energy transfer in meV
%          ydat - vector - The spectra intensity (arbitrary units)
%          edat - vector - The error in the spectra intensity (a.u.)
%          T    - scalar - The temperature at which the spectra was taken (Kelvin)
%          Ei   - scalar - The incident energy of the spectra (meV)
%   [opt]  Atot - scalar - The total area/integrated intensity of the spectra (arb. units).
%   [opt]  Ael  - scalar - The area of the elastic peak (arb. units).
%
% Outputs: spec - vector - the fitted spectra which is a sum of Gaussian peaks centred at 
%                          energies transfer corresponding to CF transitions with integrated
%                          intensities proportional to the transition matrix elements and fitted
%                          widths. The function calculates spec at values of xdat.
%          pars - vector - fit parameters - [num_peaks cen1 fwhm1 area1 ... cenN fwhmN areaN]
%          Atot - scalar - The total area/integrated intensity of the spectra (arb. units).
%          Ael  - scalar - The area of the elastic peak (arb. units).
%                          Calculated by numerical integration of spline interpolated data.

% This function uses the Levenberg-Marquardt algorithm implemented by Shrager, Jutan and
%   Muzic in the file speclsqr.m from the spec1d package to test for the lineshape that 
%   best matches the elastic peak, from choice of Gaussian, Lorentzian, and a pseudo-Voigt
%   function (superposition of Gaussian and Lorentzian) with the centres and height fixed,
%   with lfrac as the ratio of Lorentzian to Gaussian FWHM.

% Duc Le - Sat Aug 11 21:05:00 BST 2007 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 

% Some parameters
fwhm = 0.2;             % The default inital FWHM for an inelastic peak

% Parses inputs
if iscell(VHcf) && length(VHcf)==4                     % Cell array of J and CF parameters
  if isscalar(VHcf{1}) && ~mod((VHcf{1}*2),1)
    J = VHcf{1};
  else
    error('If input VHcf is cell, first cell must be scalar integer/half-integer J');
  end
  if isvector(VHcf{2}) && length(VHcf{2})==5
    V2 = VHcf{2};
  else
    error('If input VHcf is cell, second cell must be vector length 5, V2');
  end
  if isvector(VHcf{3}) && length(VHcf{3})==9
    V4 = VHcf{3};
  else
    error('If input VHcf is cell, second cell must be vector length 9, V4');
  end
  if isvector(VHcf{4}) && length(VHcf{4})==13
    V6 = VHcf{4};
  else
    error('If input VHcf is cell, second cell must be vector length 13, V6');
  end
  peaks = cflvls(cf_hmltn(J,V2,V4,V6),T,[0 1]);
elseif size(VHcf,1)==size(VHcf,2)                      % square matrix - CF Hamiltonian
  peaks = cflvls(VHcf,T,[0 1]);
elseif size(VHcf,2)==2                                 % two column vector peaks/matrix element pairs
  peaks = VHcf;
else
  error('Input VHcf must be a length 4 cell array, a square matrix or an Nx2 matrix.');
end

% Finds the elastic peak parameters
if ~exist('elas_area')
  [lineshape,elas_pars] = saficf_elaspars(xdat,ydat);
  elas_area = elas_pars(3);
end

% Calculates the area of the inelastic peaks by numerical integration of the spline interpolated spectra.
if ~exist('total_area')
  spline_spectra = @(xi) interp1(xdat,ydat,xi,'spline');
  total_area = quad(spline_spectra,min(xdat),max(xdat));
end

% Puts the peaks/matrix elements into energy transfer order
%for ip = 1:size(peaks,1)
%  peakpairs{i} = peaks(ip,:);
%end
%sort_en = peakpairs(sort(peaks(:,1)));
%for ip = 1:size(peaks,1)
%  sort_peaks(ip,:) = sort_en{ip};
%end

% Gets rid of all transitions higher than Ei
%peaks(find(peaks(:,1)>Ei),:) = [];

% Scales the transition matrix elements to the area of the inelastic spectra.
inelas_area = total_area - elas_area;
total_me = sum( peaks(find(peaks(:,1)<Ei),2) );       % Want only matrix elements of peaks < Ei
scale_factor = total_me / inelas_area;

% Finds the integrated area for each peak
pars = [1 elas_pars(1:3)']; const = [0 0 0 0];
for ip = 1:size(peaks,1)
  if peaks(ip,1)<Ei
    pars  = [pars peaks(ip,1) fwhm peaks(ip,2)/scale_factor];
    const = [const 0 1 0];
  end
end
pars(1) = (length(pars)-1)/3;

% Does the fitting!
fitpars = saficf_fitpeaks(xdat,ydat,edat,pars,const);
spec = saficf_peakfitfunc(xdat,fitpars);

%-------------------------------- Fit function - gaussians ------------------------------------ %

function out = saficf_peakfitfunc(data,pars)
% Calculates the fitting function of n gaussian peaks. 
n = pars(1); pars(1) = []; out{1} = zeros(size(data));
for i_pk = 1:3:(3*n)
  out{(i_pk-1)/3+2} = feval(@fgauss,data,pars(i_pk:(i_pk+2)));
  out{1} = out{1} + out{(i_pk-1)/3+2};
end

