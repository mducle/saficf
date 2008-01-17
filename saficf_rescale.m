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

% Defines the elastic range as 1.5*FWHM of the elastic peak for each dataset. The elastic
%   range is the range in energy transfer where SAFiCF will not look for inelastic peaks.
for i_set = 1:num_dataset
  elaspars = saficf_elaspars(xdat{i_set},ydat{i_set});
  elas_rng(i_set) = elaspars{2}(2) * 1.5;
end

% Find datasets sharing same T,Ei,freq and adds statistics in them!
indexcounter = ones(1,num_dataset); 
% Finds the indices of datasets with the same values of T, Ei and freq, and puts each unique 
%   set of datasets into a cell.
for i_set = 1:num_dataset 
  i_match_all{i_set} = find(T==T(i_set) & Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0);
  indexcounter(find(T==T(i_set) & Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0)) = 0; 
end
% Deletes the empty cells in the array i_match. 
i_match_all(find(cellfun('isempty',i_match_all))) = [];
% Loops over cells in i_match_all and adds statistics of datasets in each cell.
for i_cell = 1:length(i_match_all)
% Make sure that datasets with same T, Ei, freq have same binning to add statistics!
  n_sets_sum = length(i_match_all{i_cell});
  if sum( sum(cell2mat(xdat(i_match_all{i_cell})),2) ./ n_sets_sum ... 
            - xdat{i_match_all{i_cell}(1)} ) < 1e-3
    ydat{i_match_all{i_cell}(1)} = sum(cell2mat(ydat(i_match_all{i_cell})),2) ./ n_sets_sum;
    if ~isempty(edat)
      edat{i_match_all{i_cell}(1)} = sum(cell2mat(edat(i_match_all{i_cell})),2) ./ n_sets_sum;
    end
  end
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
    first_set_x = xdat{i_match_Ei{i_cell}(1)};
    first_set_y = ydat{i_match_Ei{i_cell}(1)};
    ref_elas_height = max( first_set_y( find(abs(first_set_x) < elas_rng(i_set)) ) );
    for i_set = i_match_Ei{i_cell}(2:length(i_match_Ei{i_cell}))
      ratio_elas_heights = ref_elas_height / max(ydat{i_set}(find(abs(xdat{i_set})<elas_rng(i_set))));
      ydat{i_set} = ydat{i_set} .* ratio_elas_heights;
      if ~isempty(edat)
        edat{i_set} = edat{i_set} .* ratio_elas_heights;
      end
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
    ref_inelas_ht = max(ydat{i_match_T{i_cell}(1)}(find(xdat{i_match_T{i_cell}(1)}>elas_rng(i_set))));
    ref_inelas_Et = xdat{i_match_T{i_cell}(1)}(find(ydat{i_match_T{i_cell}(1)}==ref_inelas_ht));
    % Checks that ref_inelas_Et is scalar!
    if ~isscalar(ref_inelas_Et)
      % Checks that the values of ref_inelas_Et are similar. peak_tol is given in saficf_defs.m
      if ( ref_inelas_Et(1) - ( sum(ref_inelas_Et)/length(ref_inelas_Et) ) ) < peak_tol
        ref_inelas_Et = ref_inelas_Et(1);
      end
    end
    %TODO: Possible bug here if there are two peaks of roughly equal intensity.
    for i_set = i_match_T{i_cell}(2:length(i_match_T{i_cell}))
      new_inelas_ht = max(ydat{i_set}(find(xdat{i_set}>elas_rng(i_set))));
      % Checks that we are using the same peak!
      if abs(ref_inelas_Et-xdat{i_set}(find(ydat{i_set}==new_inelas_ht))) < peak_tol
        ratio_inelas_ht = ref_inelas_ht / new_inelas_ht;
        ydat{i_set} = ydat{i_set} .* ratio_inelas_ht;
        if ~isempty(edat)
          edat{i_set} = edat{i_set} .* ratio_inelas_ht;
        end
      end
    end
  end
end

