function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)

% -----------------------------  Some default parameters  ------------------------------- %
maxsplit = 50;  % meV - the full CF splitting of the spin-orbit degenrate ground state.
markov_l = 100; % Length of Markov chain.
% --------------------------------------------------------------------------------------- %

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

%max_y = [];
cost = 0;
for i_set = 1:num_dataset
  % Gets the parameters of the elastic peak and lineshape for each dataset
  [lineshape{i_set},elas_pars{i_set}] = saficf_elaspars(xdat{i_set},ydat{i_set});
  % Evaluates elastic peaks.for each dataset
  elas_peak{i_set} = feval(lineshape{i_set},xdat{i_set},elas_pars{i_set}); 
  %max_y = [max_y max(ydat(:,i_set))];
  %intfac = rand(1,num_set) .* max_y ./ 10;            % Intensity factor for inelastic peaks
  intfac{i_set} = rand*max(ydat{i_set})/10;         % Intensity factor for inelastic peaks
  spectmp = saficf_genspec(J,T(i_set),V,xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
  spec{i_set} = elas_peak{i_set} + spectmp.*intfac{i_set};

  % Draws initial generated spectra onto respective axes.
  if exist('handle') && ishandle(handle(i_set))
    xi = min(xdat{i_set}):0.1:max(xdat{i_set});
    spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape{i_set});
    e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
    spcp = e_pk(:) + spcp(:).*intfac{i_set};
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

while stopflag < 5                                  % Terminates schedule after 5 cycles
  transflag = 0;                                    %   where energy(cost) has not changed
  disp([cost c intfac{1} intfac{2}]);
  for ind_markov = 1:markov_l                       % Set markov chain length arbitrarily
    cost_new = 0;                                   % Resets cost value to sum over datasets
    Vnew     = saficf_perturb(V,c,range);           % Generates new configuration
    for i_set = 1:num_dataset                       %   and new intensity factor
      intfac_n{i_set} = abs( intfac{i_set} + lrnd(c)/10 );
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
      spcp = e_pk(:) + spcp(:).*intfac{i_set};
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
