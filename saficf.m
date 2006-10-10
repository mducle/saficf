function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat)

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
  intfac{iset} = rand*max(ydat{i_set})/10;          % Intensity factor for inelastic peaks
  spectmp = saficf_genspec(J,T(i_set),V,xdat{i_set},Ei(i_set),freq(i_set),lineshape);
  spec{i_set} = elas_peak{i_set} + spectmp.*intfac;

  % Draws initial generated spectra onto respective axes.
  if exist('handle') && ishandle(handle(i_set))
    xi = min(xdat{i_set}):0.1:max(xdat{i_set});
    spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape);
    e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
    spcp = e_pk(:) + spcp(:).*intfac;
    hold all;
    plot(handle(i_set),xi,spcp,'-r');
    hold off;
    drawnow;
  end

  % Calculates the initial cost
  if max(max(abs(edat{i_set}))) < 1e-5
    cost = cost + sqrt( sum( (spec{i_set} - ydat{i_set}).^2 ) );
    costflag = 1;                                   % Cost is root mean square difference
  else
    cost = cost + sum( (spec{i_set} - ydat{i_set}).^2 ./ (edat{i_set}.^2) ) ;
    costflag = 0;                                   % Cost is chi-square
  end

end

stopflag = 0;
iteration = 1;

cstart = saficf_startc(xdat,ydat,J,T,Ei,freq,ptgpstr,elas_peak);
c = cstart;
disp('  Cost        c         intfac');

while stopflag < 5                                  % Terminates schedule after 5 cycles
  transflag = 0;                                    %   where energy(cost) has not changed
  disp([cost c intfac]);
  for ind_markov = 1:markov_l                       % Set markov chain length arbitrarily
    Vnew     = saficf_perturb(V,c,range);           % Generates new configuration
    for i_set = 1:num_set                           %   and new intensity factor
      intfac_n{i_set} = abs( intfac{i_set} + lrnd(c)/10 );
      spectmp = saficf_genspec(J,T(i_set),V, ...
        xdat{i_set},Ei(i_set),freq(i_set),lineshape);
      spec_new = elas_peak{i_set} ...               % Add elastic and inelastic peaks
                 + spectmp.*intfac_n{i_set};
      cost_new = cost;                              % Calculates new cost 
      if costflag
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
  for i_set = 1:num_set
    if exist('handle') && ishandle(handle(i_set))   % Updates graphs
      xi = min(xdat{i_set}):0.1:max(xdat{i_set});
      spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape);
      e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
      spcp = e_pk(:) + spcp(:).*intfac;
      hold all;
      plot(handle(i_set),xi,spcp,'-r');
      hold off;
      drawnow;
    end
  end
end

[cstart c iteration intfac] 

% Draws final graphs
for i_set = 1:num_set
  if exist('handle') && ishandle(handle(i_set))
    xi = min(xdat{i_set}):0.1:max(xdat{i_set});
    spcp = saficf_genspec(J,T(i_set),V,xi,Ei(i_set),freq(i_set),lineshape);
    e_pk = feval(lineshape{i_set},xi,elas_pars{i_set});
    spcp = e_pk(:) + spcp(:).*intfac;
    hold all;
    plot(handle(i_set),xi,spcp,'-r');
    hold off;
    drawnow;
  end
end
