function fitpars = saficf_fitpeaks(xdat,ydat,edat,grph,const)
% Returns the parameters [x1 fwhm1 area1 ... xn fwhmn arean] for the fit to a set of (x,y,e) data. 
%
% Inputs:  xdat  - vector - the x-data
%          ydat  - vector - the y-data
%          edat  - vector - the errors in the data
%          grph  - vector - flag (scalar logical=1 to plot graph) to specify graphical determination
%                           of peak positions or estimated parameters (vector). A scalar non logical
%                           1 means graph is assumed to be plotted already.
%          const - vector - constraints (if any) [NB. don't use if grph is not estimated parameters]
%
% Outputs: peaks - vector - [number_of_peaks x1 fwhm1 area1 ... xn fwhmn arean]

% This function uses the Levenberg-Marquardt algorithm implemented by Shrager, Jutan and
%   Muzic in the file speclsqr.m from the spec1d package to test for the lineshape that 
%   best matches the elastic peak, from choice of Gaussian, Lorentzian, and a pseudo-Voigt
%   function (superposition of Gaussian and Lorentzian) with the centres and height fixed,
%   with lfrac as the ratio of Lorentzian to Gaussian FWHM.

% Duc Le - Sat Aug  4 01:40:23 BST 2007 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 

% Checks input are all vectors
if ~isvector(xdat) | ~isvector(ydat) | ~isvector(edat)
  error('Input xdat, ydat, edat must all be vectors');
end

% Puts x,y,e into column format
xdat = xdat(:);
ydat = ydat(:);
edat = edat(:);

% Estimates starting parameters of fit ...
if exist('grph')                                        % ... by user clicking on graph.
  if isscalar(grph) && (grph==1)                        % If grph is not logical 1, assume
    errorbar(xdat,ydat,edat);                           %   user has already plotted data.
  elseif isvector(grph) && length(grph)==(grph(1)*3+1)  % If grph is vector assume it is
    estpars = grph;                                     %   the estimated parameters.
    n = estpars(1);
  else
    finish = 0; n = 0; estpars = 1;
    while (finish==0)
      n = n + 1; 
      disp(['Click on peak ' num2str(n) '. Right click to finish.'])
      waitforbuttonpress;
      [finish,peaks(n,:)] = saficf_buttonpress();
      if finish==1                                      % If use right clicks now, don't
        peaks(n,:) = []; n = n - 1;                     %   save this point and break out
        break                                           %   of the loop.
      end
      disp(['Click on half width for peak ' num2str(n)])  
      waitforbuttonpress;
      [finish,fwhmxye] = saficf_buttonpress(); 
      fwhm = abs(peaks(n,1)-fwhmxye(1)) * 2;
      estpars = [estpars peaks(n,1) fwhm fwhm*peaks(n,2)];
      estpars(1) = n;
    end
    xpeaks = peaks(:,1); ypeaks = peaks(:,2);
  end
else                                                    % ... from data directly using
  [n,xpeaks,ypeaks] = saficf_findpeaks(xdat,ydat,edat); %     _findpeaks and chopem.
  Ei = max(xdat)+5;
  estpars = n; 
  for i_pk = 1:n
    fwhm = chopem(xpeaks(i_pk),Ei,200); 
    estpars = [estpars xpeaks(i_pk) fwhm fwhm*ypeaks(i_pk)];
  end
end

if exist('const') && length(const)==length(estpars)
  warning('off','MATLAB:conversionToLogical');          % Stops matlab complaining about
  logconst = logical(const); [wnmsg,wncode] = lastwarn; %   converting numbers to 1s
  if ~strcmp(wncode,'MATLAB:conversionToLogical')       % If const is logical or 1s, times
    const = [0 ones(1,3*n)*0.1].*logconst;              %   it by 0.1 for dpin. Else assume
  end                                                   %   const==dpin.
else
  const = [0 ones(1,3*n)*0,1];
end

%[fitpars,stdev] = speclsqr(xdat,ydat,edat,estpars,const,@saficf_peakfitfunc);
[fitpars,stdev] = saficf_salsqr(xdat,ydat,edat,estpars,const,@saficf_peakfitfunc);

figure; errorbar(xdat,ydat,edat); hold all; plot(xdat,saficf_peakfitfunc(xdat,fitpars)); hold off;

%-------------------------------- Fit function - gaussians ------------------------------------ %

function out = saficf_peakfitfunc(data,pars)
% Calculates the fitting function of n gaussian peaks. 
n = pars(1); pars(1) = []; out = zeros(size(data));
for i_pk = 1:3:(3*n)
  out = out + feval(@fgauss,data,pars(i_pk:(i_pk+2)));
end

