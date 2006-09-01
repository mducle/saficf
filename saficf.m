function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat)

% -----------------------------  Some default parameters  ------------------------------- %
maxsplit = 50; % meV - the full CF splitting of the spin-orbit degenrate ground state.
% --------------------------------------------------------------------------------------- %

if ~exist('xdat')
  [xdat,ydat,edat,handle] = saficf_getdata;
end

xdat = xdat(:); ydat = ydat(:);

stopflag = 0;
iteration = 1;

elas_pars = saficf_elaspars(xdat,ydat);             % Gets parameters of elastic peak
lineshape = elas_pars{1};                           % Sets lineshape for all peaks
elas_peak = feval(lineshape,xdat,elas_pars{2});     % Evaluates elastic peaks.

intfac = rand*max(ydat)/10;                         % Intensity factor for inelastic peaks
V = saficf_genstart(ptgpstr,maxsplit,J);
spec = saficf_genspec(J,T,V,xdat,Ei,freq,lineshape);
spec = elas_peak + spec.*intfac;

if exist('handle')
  xi = min(xdat):0.1:max(xdat);
  spcp = saficf_genspec(J,T,V,xi,Ei,freq,lineshape);
  e_pk = feval(lineshape,xi,elas_pars{2});
  spcp = e_pk(:) + spcp(:).*intfac;
  hold all;
  old_plot_handle = plot(xi,spcp,'-r');
  hold off;
  drawnow;
end

cost = sqrt( sum( (spec - ydat).^2 ) );

cstart = saficf_startc(xdat,ydat,J,T,Ei,freq,ptgpstr,elas_peak);
c = cstart;
disp('  Cost        c         intfac');

while stopflag < 5                                  % Terminates schedule after 5 cycles
  transflag = 0;                                    %   where energy(cost) has not changed
  disp([cost c intfac]);
  for ind_markov = 1:100                            % Set markov chain length arbitrarily
    Vnew     = saficf_perturb(V,c);                 % Generates new configuration
    intfac_n = abs( intfac + lrnd(c) );             %   and new intensity factor
    spec_new = saficf_genspec(J,T,V,xdat,Ei,freq,...
               lineshape);
    spec_new = elas_peak + spec_new.*intfac;        % Add elastic and inelastic peaks
    cost_new = sqrt( sum( (spec_new - ydat).^2 ) ); % Cost is root mean square difference
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
  if c/cstart < 1e-3 || iteration > 1000
    stopflag = 6;
  end
  if exist('handle')
    xi = min(xdat):0.1:max(xdat);
    spcp = saficf_genspec(J,T,V,xi,Ei,freq,lineshape);
    e_pk = feval(lineshape,xi,elas_pars{2});
    spcp = e_pk(:) + spcp(:).*intfac;
    delete(old_plot_handle);
    hold all;
    old_plot_handle = plot(xi,spcp,'-r');
    hold off;
    drawnow;
  end
end

[cstart c iteration intfac] 

if exist('handle')
  xi = min(xdat):0.1:max(xdat);
  spcp = saficf_genspec(J,T,V,xi,Ei,freq,lineshape);
  e_pk = feval(lineshape,xi,elas_pars{2});
  spcp = e_pk(:) + spcp(:).*intfac;
  hold all;
  plot(xi,spcp);
  hold off;
end
