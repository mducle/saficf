function peaks = saficf_findpeaks(xdat,ydat,edat)
% Returns the x-coordinate of the peaks in a set of (x,y,e) data. 
%
% Inputs:  xdat  - vector     - the x-data
%          ydat  - vector     - the y-data
%          edat  - vector     - the errors in the data
%
% Outputs: peaks - cell array - {number_of_peaks, x-coordinates_of_peaks}

% By Duc Le 2006 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 
% Fri Oct  6 16:21:31 BST 2006 - First version.

% Puts x,y,e into column format
xdat = xdat(:);
ydat = ydat(:);
edat = edat(:);

% Finds the numerical first derivatives of the data.
dy_dx = diff(ydat) ./ diff(xdat);

% For real data need to bin first derivative to get smoother curve.
if mod(length(dy_dx),2) == 0   % Even length
  dydx_binned = ( dy_dx(1:2:length(dy_dx)) + dy_dx(2:2:length(dy_dx)) )./2;
  x_binned = ( xdat(1:2:length(xdat)-1) + xdat(2:2:length(xdat)) )./2;
  % Finds the numerical second derivative of the data.
  d2y_dx2 = diff(dydx_binned) ./ diff(x_binned);
else                           % Odd length
  dydx_binned = ( dy_dx(1:2:length(dy_dx)-1) + dy_dx(2:2:length(dy_dx)) )./2;
  x_binned = ( xdat(1:2:length(xdat)) + xdat(2:2:length(xdat)) )./2;
  % Finds the numerical second derivative of the data.
  d2y_dx2 = diff(dydx_binned) ./ diff(x_binned(2:length(x_binned)));
end

% Finds the numerical second derivative of the data.
%d2y_dx2 = diff(dy_dx) ./ diff(xdat(2:length(xdat)));
%d2y_dx2 = diff(dydx_binned) ./ diff(x_binned(2:length(x_binned)));
%d2y_dx2 = diff(dydx_binned) ./ diff(x_binned);

% Finds index of turning points where dy/dx=0
dydxzero_ind = [];
for ind_x = 2:length(x_binned)-2
  %if (dy_dx(ind_x)>0 & dy_dx(ind_x-1)<0) || (dy_dx(ind_x)<0 & dy_dx(ind_x-1)>0)
  if (dydx_binned(ind_x)>0 & dydx_binned(ind_x-1)<0) || ...
     (dydx_binned(ind_x)<0 & dydx_binned(ind_x-1)>0)
    dydxzero_ind = [dydxzero_ind ind_x];
  end
end

% Finds the magnitude of the elastic peak
elas_mag = max(ydat(find(abs(xdat)<2)));

% Finds the peak position index, from maximum in spectra where dy/dx=0 and d2y/dx2<0
peak_pos = dydxzero_ind(find(d2y_dx2(dydxzero_ind)<0));
threshold = elas_mag/50;
peak_pos = peak_pos(find(d2y_dx2(peak_pos)<-threshold));

peaks = {length(peak_pos), x_binned(peak_pos)};
