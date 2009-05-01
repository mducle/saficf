function [n,px,py] = saficf_findpeaks(xdat,ydat,edat)
% Returns the x- and y-coordinates of the peaks in a set of (x,y,e) data. 
%
% Inputs:  xdat  - vector     - the x-data
%          ydat  - vector     - the y-data
%          edat  - vector     - the errors in the data
%
% Outputs: peaks - cell array - {number_of_peaks, x-coordinates_of_peaks}

% By Duc Le 2006 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 
% Fri Oct  6 16:21:31 BST 2006 - First version.

% mdl - updated 070206: added y-coordinate output

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
  y_binned = ( ydat(1:2:length(ydat)-1) + ydat(2:2:length(ydat)) )./2;
  % Finds the numerical second derivative of the data.
  d2y_dx2 = diff(dydx_binned) ./ diff(x_binned);
else                           % Odd length
  dydx_binned = ( dy_dx(1:2:length(dy_dx)-1) + dy_dx(2:2:length(dy_dx)) )./2;
  x_binned = ( xdat(1:2:length(xdat)) + xdat(2:2:length(xdat)) )./2;
  y_binned = ( ydat(1:2:length(ydat)) + ydat(2:2:length(ydat)) )./2;
  % Finds the numerical second derivative of the data.
  d2y_dx2 = diff(dydx_binned) ./ diff(x_binned(2:length(x_binned)));
end

% Finds index of turning points where dy/dx=0
dydxzero_ind = [];
for ind_x = 2:length(x_binned)-2
  if (dydx_binned(ind_x)>0 & dydx_binned(ind_x-1)<0) || ...
     (dydx_binned(ind_x)<0 & dydx_binned(ind_x-1)>0)
    dydxzero_ind = [dydxzero_ind ind_x];
  end
end

% Finds the magnitude of the elastic peak
elas_mag = max(ydat(find(abs(xdat)<2))); if isempty(elas_mag); elas_mag = max(ydat); end;

% Finds the peak position index, from maximum in spectra where dy/dx=0 and d2y/dx2<0
peak_pos = dydxzero_ind(find(d2y_dx2(dydxzero_ind)<0));
threshold = elas_mag/50;
peak_pos = peak_pos(find(d2y_dx2(peak_pos)<-threshold));

if nargin == 1
  n = {length(peak_pos), x_binned(peak_pos) y_binned(peak_pos)};
else
  n = length(peak_pos);
  px = x_binned(peak_pos);
  py = y_binned(peak_pos);
end
