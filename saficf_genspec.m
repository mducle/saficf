function [spec,peaks] = saficf_genspec(J,T,V,Et,Ei,freq,lineshape)
% Generates an inelastic neutron spectra from a set of Fabi normalised CF parameters
%
% Syntax:  spec = saficf_genspec(J,T,V,Et,Ei,freq)
%
% Inputs:  J         - scalar     - total angular momentum quantum number
%          T         - scalar     - temperature of sample (Kelvin)
%          V         - cell array - {V2 V4 V6} CF parameters
%          Et        - vector     - energy transfer to calculate at (meV)
%          Ei        - scalar     - incident energy (meV)
%          freq      - scalar     - chopper frequency (Hz)
%          lineshape - string     - the lineshape type: 'fgauss','florr' or 'fvoigt'
%
% Outputs: spec      - vector     - the inelastic spectrum evaluated at Et points.

%-----------------------------------------------------------------------------------%
f = 0.1;   % Lorentzian factor for pseudo-voigt lineshape. 
%-----------------------------------------------------------------------------------%
%TODO: Fix so that can vary/fit lorentzian factor rather than have constant as here

peaks = [];
for ind_sites = 1:size(V,2)
  Hcf = norm_cfhmltn(J,V(:,ind_sites));
  peaks = [peaks; cflvls(Hcf,T,[0 1])];
end

fwhm = chopem(peaks(:,1),Ei,freq);

for ind_p = 1:size(peaks,1)
  peak(:,ind_p) = feval(lineshape,Et,[peaks(ind_p,1) fwhm(ind_p) peaks(ind_p,2) f]);
end

for ind_Et = 1:length(Et)
  spec(ind_Et) = sum(peak(ind_Et,:));
end

%spec = sum_peaks;

%function spec = saficf_genspec(J,T,V,Et,Ei,freq)
%
%%Hcf = cf_hmltn(J,B2,B4,B6);
%Hcf = norm_cfhmltn(J,V);
%peaks = cflvls(Hcf,T,[0 1]);
%
%%fwhm = rand(1,size(peaks,1)).*6;
%fwhm = chopem(peaks(:,1),Ei,freq);
%
%sigma = fwhm ./ (2*sqrt(2*log(2)));
%
%for ind_p = 1:size(peaks,1)
%  for ind_Et = 1:size(Et,2)
%    gauss_peak(ind_Et,ind_p) = peaks(ind_p,2) * exp(-(Et(ind_Et)-peaks(ind_p,1))^2 ...
%                                                     ./ (2.*sigma(ind_p)^2));
%  end
%end
%
%for ind_Et = 1:size(Et,2)
%  sum_peaks(ind_Et) = sum(gauss_peak(ind_Et,:));
%end
%
%%spectra = [Et' sum_peaks'];
%%spec = Et' + i.*sum_peaks';
%spec = sum_peaks;
