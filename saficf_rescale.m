function [T,Ei,freq,xdat,ydat,edat] = saficf_rescale(T,Ei,freq,xdat,ydat,edat)
% saficf_rescale - rescales data sets with same temperature or Ei, chopper freq.
%
% Syntax:  [T,Ei,freq,xdat,ydat,edat] = saficf_rescale(T,Ei,freq,xdat,ydat,edat)
%
% Inputs:  T    - vector     - the temperatures of the datasets
%          Ei   - vector     - the incident energy of the datasets
%          freq - vector     - the chopper frequency of the datasets
%          xdat - cell array - the energy transfers as a vector in each cell 
%          ydat - cell array - the intensity/scattering function as a vector per cell. 
%          edat - cell array - the error in the y-data as a vector per cell.
%
% Outputs: xdat - cell array - the energy transfer, copied from the input
%          ydat - cell array - the intensity datasets rescaled if same T, Ei/freq
%          edat - cell array - the error in y, rescaled if same T, Ei/freq.

% This file is part of the SAFiCF package. Licensed under the Gnu GPL v2.

% Duc Le - duc.le@ucl.ac.uk - Mon Oct 23 11:50:24 BST 2006 

% Include definitions
saficf_defs;

num_dataset = length(T);

% Find datasets sharing same T,Ei,freq and adds statistics in them!
indexcounter = ones(1,num_dataset); 
% Finds the indices of datasets with the same values of T, Ei and freq, and puts each unique 
%   set of datasets into a cell.
for i_set = 1:num_dataset; 
  i_match_all{i_set} = find(T==T(i_set) & Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0);
  indexcounter(find(T==T(i_set) & Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0)) = 0; 
end
% Deletes the empty cells in the array i_match. 
i_match_all(find(cellfun('isempty',i_match_all))) = [];
% Loops over cells in i_match_all and adds statistics of datasets in each cell.
for i_cell = 1:length(i_match_all)
%TODO: make sure that datasets with same T, Ei, freq have same binning to add statistics!
  n_sets_sum = length(i_match_all{i_cell});
  ydat{i_match_all{i_cell}(1)} = sum(cell2mat(ydat(i_match_all{i_cell})),2) ./ n_sets_sum;
  edat{i_match_all{i_cell}(1)} = sum(cell2mat(edat(i_match_all{i_cell})),2) ./ n_sets_sum;
end
% Deletes other cells that have now been summed
i_match_all = cell2mat(i_match_all'); 
now_empty = i_match_all(:,2:n_sets_sum)';
xdat(now_empty) = [];
ydat(now_empty) = [];
edat(now_empty) = [];
T(now_empty) = [];
Ei(now_empty) = [];
freq(now_empty) = [];
num_dataset = length(T);

% Scales the data in datasets with same Ei,freq but different T to the elastic peak height
indexcounter = ones(1,num_dataset); 
% Finds the indices of datasets with the same values of Ei and freq, and puts each unique 
%   set of datasets into a cell.
for i_set = 1:num_dataset; 
  i_match_Ei{i_set} = find(Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0);
  indexcounter(find(Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0)) = 0; 
end
% Deletes the empty cells in the array i_match. 
i_match_Ei(find(cellfun('isempty',i_match_Ei))) = [];
% Loops over cells in i_match_Ei and scales datasets with respects to first dataset in cell.
for i_cell = 1:length(i_match_Ei)
  if length(i_match_Ei{i_cell}) ~= 1   % Leaves datasets with unique T,Ei,freq alone. 
    ref_elas_height = max(ydat{i_match_Ei{i_cell}(1)}(find(abs(xdat{i_match_Ei{i_cell}(1)})<elas_rng)));
    for i_set = i_match_Ei{i_cell}(2:length(i_match_Ei{i_cell}))
      ratio_elas_heights = ref_elas_height / max(ydat{i_set}(find(abs(xdat{i_set})<elas_rng)));
      ydat{i_set} = ydat{i_set} .* ratio_elas_heights;
      edat{i_set} = edat{i_set} .* ratio_elas_heights;
    end
  end
end

% Scales same-temperature, different Ei or freq datasets to largest inelastic peak height.
indexcounter = ones(1,num_dataset); 
% Finds indices of datasets with same value of T and put each set of datasets into a cell.
for i_set = 1:num_dataset; 
  i_match_T{i_set} = find(T==T(i_set) & indexcounter~=0);
  indexcounter(find(T==T(i_set) & indexcounter~=0)) = 0; 
end
% Deletes the empty cells in the array i_match. 
i_match_T(find(cellfun('isempty',i_match_T))) = [];
% Loops over cells in i_match_T and scales datasets with respect to largest inelastic peak.
for i_cell = 1:length(i_match_T)
  if length(i_match_T{i_cell}) ~= 1    % Leaves datasets with unique T,Ei,freq alone. 
    ref_inelas_ht = max(ydat{i_match_T{i_cell}(1)}(find(xdat{i_match_T{i_cell}(1)}>elas_rng)));
    ref_inelas_Et = xdat{i_match_T{i_cell}(1)}(find(ydat{i_match_T{i_cell}(1)}==ref_inelas_ht));
      %TODO: Check that ref_inelas_Et is scalar!
    for i_set = i_match_T{i_cell}(2:length(i_match_T{i_cell}))
      new_inelas_ht = max(ydat{i_set}(find(xdat{i_set}>elas_rng)));
      % Checks that we are using the same peak!
      if abs(ref_inelas_Et-xdat{i_set}(find(ydat{i_set}==new_inelas_ht))) < 0.1
        ratio_inelas_ht = ref_inelas_ht / new_inelas_ht;
        ydat{i_set} = ydat{i_set} .* ratio_inelas_ht;
        edat{i_set} = edat{i_set} .* ratio_inelas_ht;
      end
    end
  end
end
%TODO: Change criterion so that elas_rng = 1.5*FWHM

