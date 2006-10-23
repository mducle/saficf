function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)

% -----------------------------  Some default parameters  ------------------------------- %
maxsplit = 50;  % meV - the full CF splitting of the spin-orbit degenrate ground state.
elas_rng = 2;   % meV - the range from +elas_rng to -elas_rng in Et for the elastic peak.
markov_l = 100; % Length of Markov chain.
% --------------------------------------------------------------------------------------- %

%TODO: Put all input checking into another file, e.g. saficf_inputs.m
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
if ~exist('xdat')
  for i_set = 1:num_dataset
    [xdat{i_set},ydat{i_set},edat{i_set},handle(i_set)] = saficf_getdata;
    xdat{i_set} = xdat{i_set}(:); 
    ydat{i_set} = ydat{i_set}(:);
    edat{i_set} = edat{i_set}(:);
  end
else
  if ~exist('ydat')
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
    if exist('edat')
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

% Generate starting values for the CF parameters, and their max/min ranges.
V = saficf_genstart(ptgpstr,maxsplit,J);
range = saficf_range(ptgpstr,maxsplit,J);

% Sets an intensity factor for the inelastic peaks.
intfac = rand*max(ydat{1})/10;

%TODO: Put all this rescaling data into another file, e.g. saficf_rescale.m
% Scales the data in datasets with same Ei,freq but different T to the elastic peak height
indexcounter = ones(1,num_dataset); 
% Finds the indices of datasets with the same values of Ei and freq, and puts each unique 
%   set of datasets into a cell.
for i_set = 1:num_dataset; 
  i_match_Ei{i_set} = find(Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0);
  indexcounter(find(Ei==Ei(i_set) & freq==freq(i_set) & indexcounter~=0)) = 0; 
end
% Deletes the empty cells in the array i_match. 
i_match_Ei(find(cellfun('isempty',i_match))) = [];
% Loops over cells in i_match_Ei and scales datasets with respects to first dataset in cell.
for i_cell = 1:length(i_match_Ei)
  if length(i_match_Ei{i_cell}) ~= 1   % Leaves datasets with unique T,Ei,freq alone. 
    ref_elas_height = max(ydat{i_match_Ei{i_cell}(1)}(find(abs(xdat{i_match_Ei{i_cell}(1))<elas_rng)));
    for i_set = i_match_Ei{i_cell}(2:length(i_match_Ei{i_cell}))
      ratio_elas_heights = ref_elas_height / max(ydat{i_set}(find(abs(xdat{i_set})<elas_rng));
      ydat{i_set} = ydat{i_set} .* ratio_elas_heights;
      edat{i_set} = edat{i_set} .* ratio_elas_heights;
    end
  end
end
% TODO: Find datasets sharing same T,Ei,freq and adds statistics in them!

% Scales same-temperature, different Ei or freq datasets to largest inelastic peak height.
indexcounter = ones(1,num_dataset); 
% Finds indices of datasets with same value of T and put each set of datasets into a cell.
for i_set = 1:num_dataset; 
  i_match_T{i_set} = find(T==T(i_set) & indexcounter~=0);
  indexcounter(find(T==T(i_set) & indexcounter~=0)) = 0; 
end
% Loops over cells in i_match_T and scales datasets with respect to largest inelastic peak.
for i_cell = 1:length(i_match_T)
  if length(i_match_T{i_cell}) ~= 1    % Leaves datasets with unique T,Ei,freq alone. 
    ref_inelas_ht = max(ydat{i_match_T{i_cell}(1)}(find(xdat{i_match_T{i_cell}(1)>elas_rng));
    ref_inelas_Et = xdat{i_match_T{i_cell}(1)}(find(ydat{i_match_T{i_cell}(1)}==ref_inelas_ht));
      %TODO: Check that ref_inelas_Et is scalar!
    for i_set = i_match_T{i_cell}(2:length(i_match_T{i_cell}))
      new_inelas_ht = max(ydat{i_set}(find(xdat{i_set}>elas_rng));
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

cost = 0;
for i_set = 1:num_dataset
  % Gets the parameters of the elastic peak and lineshape for each dataset
  [lineshape{i_set},elas_pars{i_set}] = saficf_elaspars(xdat{i_set},ydat{i_set});
  % Evaluates elastic peaks.for each dataset
  elas_peak{i_set} = feval(lineshape{i_set},xdat{i_set},elas_pars{i_set}); 
  spectmp = saficf_genspec(J,T(i_set),V,xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
  spec{i_set} = elas_peak{i_set} + spectmp.*intfac;

  % Draws initial generated spectra onto respective axes.
  if exist('handle') && ishandle(handle(i_set))
    xi = min(xdat{i_set}):0.1:max(xdat{i_set});
    spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape{i_set});
    e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
    spcp = e_pk(:) + spcp(:).*intfac;
    hold (handle(i_set),'all');
    plot(handle(i_set),xi,spcp,'-r');
    hold (handle(i_set),'off');
    drawnow;
  end

  % Calculates the initial cost
  if max(max(abs(edat{i_set}))) < 1e-5
    cost = cost + sqrt( sum( (spec{i_set} - ydat{i_set}).^2 ) );
    costflag(i_set) = 1;                            % Cost is root mean square difference
  else
    cost = cost + sum( (spec{i_set} - ydat{i_set}).^2 ./ (edat{i_set}.^2) ) ;
    costflag(i_set) = 0;                            % Cost is chi-square
  end

end

stopflag = 0;
iteration = 1;

cstart = saficf_startc(xdat,ydat,edat,J,T,Ei,freq,ptgpstr,elas_peak);
c = cstart;
disp('  Cost        c         intfac');

%TODO: Find a way to make sure cost, temperature c, and CF parameters are commensurable with
%      each other -> cost should be of order unity, CF pars steps of order 0.1 and c order 10.

while stopflag < 5                                  % Terminates schedule after 5 cycles
  transflag = 0;                                    %   where energy(cost) has not changed
  disp([cost c intfac{1} intfac{2}]);
  for ind_markov = 1:markov_l                       % Set markov chain length arbitrarily
    cost_new = 0;                                   % Resets cost value to sum over datasets
    Vnew     = saficf_perturb(V,c,range);           % Generates new configuration
    for i_set = 1:num_dataset                       %   and new intensity factor
      intfac_n{i_set} = abs( intfac + lrnd(c)/10 );
      spectmp = saficf_genspec(J,T(i_set),Vnew, ...
        xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
      spec_new = elas_peak{i_set} ...               % Add elastic and inelastic peaks
                 + spectmp.*intfac_n{i_set};
      if costflag(i_set)                            % Calculates new cost 
        cost_new = cost_new + sqrt( sum( (spec{i_set} - ydat{i_set}).^2 ) );
      else
        cost_new = cost_new + sum( (spec{i_set} - ydat{i_set}).^2 ./ (edat{i_set}.^2) );
      end
    end
    if saficf_accept(cost_new-cost,c)               % Acceptance criteria is Fermi-Dirac
      V = Vnew;
      cost = cost_new;
      intfac = intfac_n;
      transflag = 1;
    end
  end
  c = saficf_schedul(cstart,iteration);
  iteration = iteration + 1;
  if transflag
    stopflag = 0;
  else
    stopflag = stopflag + 1
  end
  if c/cstart < 1e-3 || iteration > 1000            % Terminates after temperature falls
    stopflag = 6;                                   %   by 1000 or 1000 iterations.
  end
  for i_set = 1:num_dataset
    if exist('handle') && ishandle(handle(i_set))   % Updates graphs
      xi = min(xdat{i_set}):0.1:max(xdat{i_set});
      spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape{i_set});
      e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
      spcp = e_pk(:) + spcp(:).*intfac;
      hold (handle(i_set),'all');
      hold all;
      plot(handle(i_set),xi,spcp,'-r');
      hold (handle(i_set),'off');
      hold off;
      drawnow;
    end
  end
end

[cstart c iteration intfac] 

% Draws final graphs
for i_set = 1:num_dataset
  if exist('handle') && ishandle(handle(i_set))
    xi = min(xdat{i_set}):0.1:max(xdat{i_set});
    spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape{i_set});
    e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
    spcp = e_pk(:) + spcp(:).*intfac;
    hold (handle(i_set),'all');
    plot(handle(i_set),xi,spcp,'-r');
    hold (handle(i_set),'off');
    drawnow;
  end
end
