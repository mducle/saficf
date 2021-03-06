function spectra = cfins(peaks,Ei,range,fwhm)
% cfins - Generates inelastic neutron spectra from CF data
% 
% Syntax: spectra = cfins(peaks,Ei,range,fwhm)
%
% Inputs:  peaks - two column matrix with the energy and intensity of the peaks
%          Ei    - scalar with the incident energy. 
%          range - vector of form: [first_energy step last_energy],
%                  or range of values at which to calculate spectra (not three values!)
%          fwhm  - scalar with the full width at half maximum of the peaks
%
% Outputs: spectra is a two column matrix with the energy and intensity of the spectra
%
% The input peaks can be generated by the function cflvls. Type 'help cflvls'
% The inputs Ei, range, fwhm and A are optional.
% 
% If invoked as:    spectra = cfins(peaks);
% The function will output the spectra from -10:0.1:30meV with a fwhm of 1meV for all 
% peaks. 
%
% If invoked as:    spectra = cfins(peaks,Ei);
% The function will output the spectra from (-Ei/2:0.1:Ei) with fwhm determined by the 
% chopem function. 
%
% If invoked as:    spectra = cfins(peaks,[0 1 20]); or
%                   spectra = cfins(peaks,0:1:20);
% The function will output the spectra at the specified range (0:1:20) with fwhm of
% 1meV for all peaks. 
%
% If invoked as:    spectra = cfins(peaks,[0 1 20],2)
% The function will output the spectra from 0:1:20meV with fwhm of 2meV.

% Duc Le 2006 - Thu Aug 10 00:33:35 BST 2006 - duc.le@ucl.ac.uk
% mdl - updated 070126: Changed to using fgauss.m for peaks, cleaned up help section.

% This file is part of the SAfiCF package. Licenced under the GNU GPL v2 or later. 

% Checks arguments
if ~exist('peaks')
  error('You must supply a two column matrix of the energy and intensity of peaks');
end

% Cleans up peaks data from what cflvls generates
peaks = peaks(find(abs(peaks(:,2))>1e-3),:);

% Case 1: spectra = cfins(peaks);
if ~exist('Ei') 
  Et = -10:0.1:30;
  fwhm = 1;
% Case 2a: spectra = cfins(peaks,Ei);
elseif ~exist('range','var') && isscalar(Ei)
  Et = -Ei/2:0.1:Ei;
  fwhm = chopem(peaks(:,1),Ei,200);
% Case 2b: spectra = cfins(peaks,Ei,fwhm);
elseif ~exist('fwhm','var') && isscalar(Ei) && isscalar(range)
  Et = -Ei/2:0.1:Ei;
  fwhm = chopem(peaks(:,1),Ei,range);
% Case 3: spectra = cfins(peaks,range);
elseif ~exist('range','var') && isvector(Ei)
  if max(size(Ei)) == 3
    Et = Ei(1):Ei(2):Ei(3);
  else
    Et = Ei;
  end
  fwhm = 1;
% Case 4: spectra = cfins(peaks,range,fwhm);
elseif exist('range','var') && isscalar(range)
  if max(size(Ei)) == 3
    Et = Ei(1):Ei(2):Ei(3);
  else
    Et = Ei;
  end
  fwhm = range;
else
  error(['Arguments not recognised. Type ''help ' mfilename ''' for more help']);
end

% If FWHM is not determined by chopem
if isscalar(fwhm)
% Converts from FWHM to standard deviation.
  sigma = fwhm / (2*sqrt(2*log(2)));

  for ind_p = 1:size(peaks,1)
    for ind_Et = 1:size(Et,2)
      gauss_peak(:,ind_p) = feval('fgauss',Et,[peaks(ind_p,1) fwhm peaks(ind_p,2)]);
    end
  end

% FWHM is determined by chopem
else
  sigma = fwhm ./ (2*sqrt(2*log(2)));
  
  for ind_p = 1:size(peaks,1)
    for ind_Et = 1:size(Et,2)
      gauss_peak(:,ind_p) = feval('fgauss',Et,[peaks(ind_p,1) fwhm(ind_p) peaks(ind_p,2)]);
    end
  end

end % isscalar(fwhm) 

for ind_Et = 1:size(Et,2)
  sum_peaks(ind_Et) = sum(gauss_peak(ind_Et,:));
end

spectra = [Et' sum_peaks'];
