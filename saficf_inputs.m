function [xdat,ydat,edat] = saficf_inputs(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)
% Checks that the inputs to saficf are of the correct type and dimensions.
%
% Syntax:  [check,xdat,ydat,edat] = saficf_inputs(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)
%
% Inputs:  J       - scalar     - The total angular quantum number. 
%          T       - vector     - the temperature of each dataset taken.
%          Ei      - vector     - the incident energy of each dataset taken.
%          freq    - vector     - the chopper frequency of each dataset taken.
%          ptgpstr - string     - The space group symmetry string.
%          xdat    - cell array - the energy transfer as a vector, one cell per dataset.
%          ydat    - cell array - the intensity/scattering function as a vector per cell.
%          edat    - cell array - the error in y-data as a vector, one cell per dataset.
%
% Outputs: xdat    - cell array - the energy transfer datasets, from either input or from
%                                 current graphs.
%          ydat    - cell array - the intensity/scattering function datasets
%          edat    - cell array - the error in y datasets.

% This file is part of the SAFiCF package, licenced under the Gnu GPL v2.

% Sun Oct 15 20:01:06 GMT 2006 - split this section off from main saficf.m file.

% Checks input is in correct form
if ~isequal(size(J),[1 1])
  error('J must be a scalar');
elseif min(size(T)) ~= 1
  error('T must be a vector');
elseif min(size(Ei)) ~= 1
  error('Ei must be a vector');
elseif min(size(freq)) ~= 1
  error('freq must be a vector');
end  

num_dataset = length(T);
if num_dataset ~= length(Ei) || num_dataset ~= length(freq)
  error('T, Ei, freq must have the same length - i.e the number of datasets to fit');
end

% Gets data from a graph if no data provided.
if isempty(xdat)
  for i_set = 1:num_dataset
    [xdat{i_set},ydat{i_set},edat{i_set},handle(i_set)] = saficf_getdata;
    xdat{i_set} = xdat{i_set}(:); 
    ydat{i_set} = ydat{i_set}(:);
    edat{i_set} = edat{i_set}(:);
  end
else
  if isempty(ydat)
    error('If you supply x-data you must also supply the y-data!');
  elseif isvector(ydat) & isvector(xdat)
  % Convert data vectors to cell-arrays so the dataset index looping still works
    ydat = {ydat(:)};
    xdat = {xdat(:)};
    if exist('edat') 
      if isvector(edat)
        edat = {edat(:)};
      else
        error('xdat and ydat are vectors so edat must also be vector!');
      end
    end
  elseif iscell(ydat) & iscell(xdat)
    for i_set = 1:num_dataset
      xdat{i_set} = xdat{i_set}(:); 
      ydat{i_set} = ydat{i_set}(:);
    end
    if ~isempty(edat)
      if iscell(edat)
        for i_set = 1:num_dataset
          edat{i_set} = edat{i_set}(:);
        end
      else
        error('xdat and ydat are cell arrays so edat must also be a cell-array');
      end
    end
  else
    error('xdat,ydat,edat must be either vectors or cell-arrays'); 
  end 
end

