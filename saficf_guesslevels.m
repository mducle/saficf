function [energies,levels] = saficf_guesslevels(J,ptgpstr,peaks)
% Tries to assign degenerate levels for a particular point group symmetry to observed peaks in spectra
%
% Syntax:  [energies,levels] = saficf_guesslevels(J,ptgpstr,peaks)
%
% Inputs:  J        - scalar - the total angular momentum quantum number
%          ptgpstr  - string - the point group as a string (e.g. 'C2v')
%          peaks    - vector - the energy transfer of the observed peaks
%
% Outputs: energies - vector - the 2J+1 energy levels estimated
%          levels   - cell   - {[energies of singlets] [energies of doublets] ...}

% Duc Le - Tue Oct 14 15:43:43 BST 2008 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 

% Parameter that if levels are split by less than, then judge those levels are degenerate
small = 1e-5; 

% Finds the degeneracies expected for a particular point group symmetry
degen = [0 0 0 0]; levels = {[] [] [] []};       % Only go up to quartet
Hcf = cf_hmltn(J,ptgp(ptgpstr));
[V,E] = eig(Hcf);
eptgp = diag(E); eptgp = eptgp-min(eptgp); 
iL = 1;
while(eptgp)
  idegen = find(abs(eptgp-eptgp(1))<small);
  eptgp(idegen) = [];
  degen(length(idegen)) = degen(length(idegen))+1;
  lvls(iL) = length(idegen); iL=iL+1;
end

% Estimate transition matrix elements
Jmat = mag_op_j(J); Jx = Jmat(:,:,1); Jy = Jmat(:,:,2); Jz = Jmat(:,:,3);
Trans = ( (V'*Jx*V).*conj(V'*Jx*V) + (V'*Jy*V).*conj(V'*Jy*V) + (V'*Jz*V).*conj(V'*Jz*V) );

% Calculate which degenerate levels to assign to each peak
numlvls = sum(lvls);                             % Total number of levels
difflvl = numlvls-length(peaks);                 % Checks if we have unknown peaks or too many peaks!
if difflvl==0
  if(length(find(degen))==1)                     % We have all single type degeneracy 
    multiplicity = find(degen);                  %   - e.g. a Kramers ion (all doublets)
    for iE = 1:length(peaks)
      energies((iE*multiplicity):(iE*(multiplicity+1)-1)) = peaks(iE);
      levels{multiplicity} = [levels{multiplicity} peaks(iE)];
    end
  else

  end
elseif difflvl>0

elseif difflvl<0

end
