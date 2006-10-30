function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)

% Include definitions
saficf_defs;

% Checks input is in correct form
if ~exist('xdat'); xdat = {}; end
if ~exist('ydat'); ydat = {}; end
if ~exist('edat'); edat = {}; end
[xdat,ydat,edat] = saficf_inputs(J,T,Ei,freq,ptgpstr,xdat,ydat,edat);

% Rescales datasets with the same T to same elastic peak, and same Ei/freq to same inelastic peaks.
[T,Ei,freq,xdat,ydat,edat] = saficf_rescale(T,Ei,freq,xdat,ydat,edat);

num_dataset = length(T);

% Generate starting values for the CF parameters, and their max/min ranges.
V = saficf_genstart(ptgpstr,maxsplit,J);
range = saficf_range(ptgpstr,maxsplit,J);

% Sets an intensity factor for the inelastic peaks.
intfac = rand*max(ydat{1})/10;

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
  disp([cost c intfac]);
  for ind_markov = 1:markov_l                       % Set markov chain length arbitrarily
    cost_new = 0;                                   % Resets cost value to sum over datasets
    Vnew     = saficf_perturb(V,c,range);           % Generates new configuration
    intfac_n = abs( intfac + lrnd(c)/10 );          %   and new intensity factor
    for i_set = 1:num_dataset
      spectmp = saficf_genspec(J,T(i_set),Vnew, ...
        xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
      spec_new = elas_peak{i_set} ...               % Add elastic and inelastic peaks
                 + spectmp.*intfac_n;
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
