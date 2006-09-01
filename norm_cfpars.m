function [V2,V4,V6] = norm_cfpars(J, B2, B4, B6)
% norm_cfpars - generates normalised CF parameters after the method of Peter Fabi (FOCUS)
%
% Syntax:  norpars    = norm_cfpars(J, B2, B4, B6)
%     OR:  [V2,V4,V6] = norm_cfpars(J, B2, B4, B6)
%     OR:  norpars    = norm_cfpars(J, {B2 B4 B6})
%     OR:  [V2,V4,V6] = norm_cfpars(J, {B2 B4 B6})
%
% Inputs:  J   - The total angular momentum quantum number of the ground spin-orbit multiplet.
%          B2  - A 5-component vector with the crystal field parameters of rank 2. 
%          B4  - A 7-component vector with the crystal field parameters of rank 4. 
%          B6  - A 13-component vector with the crystal field parameters of rank 6. 
%          B   - A cell array {B2 B4 B6} of vectors of the crystal field parameters.
%
% Outputs: norpars - either a cell array {V2 V4 V6} containing the normalised parameters
%                    or (if invoked as the second option) three vectors with these parameters. 

% The normalisation is given by:
%                                                         2
%  q      q    q               q  2     ---    |<i|Okq|j>|
% V   =  B  ||O ||   where  ||O ||   =  >      ------------
%  k      k    k               k        ---i,j     2J+1
%
% Reference: FOCUS manual by Peter Fabi. Rutherford Appleton Laboratory Report RAL-TR-95-023 

% By Duc Le - Tue Aug 15 20:21:28 BST 2006 - duc.le@ucl.ac.uk

% This file is part of the SAfiCF package. 
% Licenced under the GNU GPL v2 or later. 

% Checks we have the correct inputs
if nargin == 2
  if ~iscell(B2)
    error(['If using only two arguments, B must be a cell array. Type "help ' mfilename '".'])
  end
elseif nargin == 4 
  if ~isvector(B2) || length(B2) ~= 5
    error(['Invalid crystal field parameters B2. Type "help ' mfilename '".'])
  elseif ~isvector(B4) || length(B4) ~= 9
    error(['Invalid crystal field parameters B4. Type "help ' mfilename '".'])
  elseif ~isvector(B6) || length(B6) ~= 13
    error(['Invalid crystal field parameters B6. Type "help ' mfilename '".'])
  end 
else
  error(['This function requires either four or two arguments. Type "help ' mfilename '".']) 
end 
if ~isscalar(J) || mod((J*2),1)
  error(['Total angular momentum J must be a scalar and integer or half integer.' ...
         ' Type "help "' mfilename '".'])
end 

% Put the parameter vectors into a cell array, so we don't have to do three loops
if iscell(B2)
  B = B2;
  [B2,B4,B6] = deal(B{:});
else
  B = {B2 B4 B6};
end

V = {zeros(1,5) zeros(1,9) zeros(1,13)};

for ind_B = 1:3
  % Don't need to convert components which are zero. So index only non zero components
  for ind_Bk = find(B{ind_B})
    Bpar_k = zeros(1,length(B{ind_B}));            Bpar_k(ind_Bk) = 1;
    Bpar   = {zeros(1,5) zeros(1,9) zeros(1,13)};  Bpar{ind_B} = Bpar_k;
    B_op = cf_hmltn(J,Bpar);
    % Long way: for i=1:(2*J+1);for j=1:(2*J+1); norm=norm+B_op(i,j)*conj(B_op(i,j));end;end
    % cf. FOCUS source, file CF_FABI.FOR line 3300
    norm = sqrt( trace(B_op*conj(B_op)) / (2*J+1) );
    V{ind_B}(ind_Bk) = B{ind_B}(ind_Bk) * norm;
  end
end

% Outputs the normalised parameters in a cell array or as seperate vectors.
if nargout == 3
  V2 = V{1}; V4 = V{2}; V6 = V{3};
else
  V2 = V;
end
