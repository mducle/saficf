function V = saficf(J,T,Ei,freq,ptgpstr,xdat,ydat,edat)

% ----------------------------------  Some default parameters  ------------------------------------ %
maxsplit = 50;  % meV - the full CF splitting of the spin-orbit degenrate ground state.
elas_rng = 3;   % meV - the range from +elas_rng to -elas_rng in Et for the elastic peak.
peak_tol = 0.1; % meV - the tolerance on the energy transfer of inelastic peak
infac_ll = 0.1; %     - the intensity factor lower limit for inelastic peaks
infac_hl = 1.5; %     - the intensity factor upper limit for inelastic peaks
maxTstep = 100; %     - maximum number of temperature steps to take
% ------------------------------------------------------------------------------------------------- %

% ------------------------------------------------------------------------------------------------- %
% Input Processing                                                                                  %
% ------------------------------------------------------------------------------------------------- %
% Checks input is in correct form
if ~exist('xdat'); xdat = {}; end
if ~exist('ydat'); ydat = {}; end
if ~exist('edat'); edat = {}; end
[xdat,ydat,edat,handle] = saficf_inputs(J,T,Ei,freq,ptgpstr,xdat,ydat,edat);

% Rescales datasets with the same T to same elastic peak, and same Ei/freq to same inelastic peaks.
[T,Ei,freq,xdat,ydat,edat] = saficf_rescale(T,Ei,freq,xdat,ydat,edat);

num_dataset = length(T);

% ------------------------------------------------------------------------------------------------- %
% Initialization                                                                                    %
% ------------------------------------------------------------------------------------------------- %

% Resets the state of the random number generator
rand('state',sum(100*clock));

% Sets an intensity factor for the inelastic peaks.
[intfac,if0] = deal(rand*max(ydat{1})/10);

% Generate starting values for the CF parameters, and their max/min ranges.
V = saficf_genstart(ptgpstr,maxsplit,J);
range = saficf_range(ptgpstr,maxsplit,J);

if (strcmp(ptgpstr,'cubic') | strncmp(ptgpstr,'o',1) | strncmp(ptgpstr,'t',1))
  n = 3;
  [x,x0] = deal([(rand(1,2)-0.5).*(2*range) intfac]); % Starting point (set of parameters)
else
  x0 = [];
  for i_s = 1:size(V,2)                               % i_s indexes sites
    %[vd(1:5),vd(6:14),vd(15:27)] = deal(V{:,i_site}); % vd is dummy variable
    vd(i_s,:) = [V{1,i_s} V{2,i_s} V{3,i_s}];         % vd is dummy variable
    nV(i_s) = length(find(vd(i_s,:)));                % Number of parameters needed per site.
    x0 = [x0 vd(i_s,find(vd(i_s,:)))];                % Starting point (set of parameters)
  end
  [x,x0] = deal([x0 intfac]); 
  n = sum(nV) + 1;                                    % Number of parameters needed.
end

% Generate a starting step vector
[v,v0] = deal( [(rand(1,n-1)-.5) rand*5].*2e-1.*x0);  % some fraction of starting point values!

% Some Simulated Annealing Parameters, after Corana et. al. 
Ns = 20;             % Number of passes to ensure step sizes give acceptance rate of approx 50%
%NT = max([100,5*n]);% Number of steps at each temperature to ensure thermal equilibrium.
NT = min([20,5*n]);  % Number of steps at each temperature to ensure thermal equilibrium.
Nepsilon = 4;        % Number of successive temperature reductions before testing for stopping fit.
ci = ones(1,n).*2;   % Step varying criterion intial value
rT = 0.85;           % Temperature reduction coefficient
epsilon = 1e-4;      % Convergence criterion

% Calculates the parameters of the elastic peak, which will not be fitted by SA, for each dataset.
cost = 0;
for i_set = 1:num_dataset
  % Gets the parameters of the elastic peak and lineshape for each dataset
  [lineshape{i_set},elas_pars{i_set}] = saficf_elaspars(xdat{i_set},ydat{i_set});
  % Evaluates elastic peaks.for each dataset
  elas_peak{i_set} = feval(lineshape{i_set},xdat{i_set},elas_pars{i_set}); 
  % Generates a spectrum with inelastic peaks from the starting CF parameters.  
  spectmp = saficf_genspec(J,T(i_set),V,xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
  spec{i_set} = elas_peak{i_set} + spectmp(:).*intfac;
  % Finds the peak positions and calculates the cost of the difference of the guess CF transition energies.
  [npeaks(i_set),xpeaks{i_set},ypeaks{i_set}] = saficf_findpeaks(xdat{i_set},ydat{i_set},edat{i_set});
  pks = [];
  for ind_sites = 1:size(V,2)
    pks = [pks; cflvls(norm_cfhmltn(4,V(:,ind_sites)),10,[0 1])]; 
  end
  pks(find(pks(:,2)<1e-2),:) = []; 
  xpks = pks(:,1); omt = [];
  for i = 1:length(xpks); 
    omt(:,i) = xpks(i)-xpeaks{i_set}; 
  end; 
  costp = sqrt(sum(min(omt.^2)));
  % Calculates the initial cost
  if isempty(edat) || max(max(abs(edat{i_set}))) < 1e-5
    cost = cost + sqrt( sum( (spec{i_set} - ydat{i_set}).^2 ) ) + costp;
    costflag(i_set) = 1;                            % Cost is root mean square difference
  else
    cost = cost + sqrt( sum( (spec{i_set} - ydat{i_set}).^2 ./ (edat{i_set}.^2) ) ) + costp;
    costflag(i_set) = 0;                            % Cost is sqrt(chi-square)
  end
  if ~isempty(handle) && ishandle(handle(i_set))
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
cost_opt = cost

% Calculates a starting temperature where at least half of all transitions are accepted.
cfpar_steps = sum(range)/length(range)/NT;          % Average change in CF par / step.
[T0,Tsa] = deal( saficf_startc(xdat,ydat,edat,J,T,Ei,freq,ptgpstr,elas_peak,cfpar_steps)*10 );

% More Simulated Annealing Parameters, after Corana et. al. 
Vopt = V; cost_opt = cost;           % Optimized CF parameters and cost function
n_u = zeros(1,n);                    % u=1,...,n; n_u is a step variation count.
f_ustar = ones(1,Nepsilon) .* cost;  % u=0,-1,...,-Nepsilon+1; stores cost of each T step.

%while ~terminating_criterion 
% ------------------------------------------------------------------------------------------------- %
% Steps 1-4                                                                                            %
% ------------------------------------------------------------------------------------------------- %
for k = 0:maxTstep
  tic
  for m = 1:NT                                             % One m-loop takes ~25min!
   %m
   %tic
   for j = 1:Ns
    for h = 1:n
      eh = zeros(1,n); eh(h) = 1;
      xprime = x + ( (rand-.5)*2*v(h) .* eh );             % Step 1
      if h<n
        if xprime(h) < -range | xprime(h) > range
          break;                                           % Step 2
        end
      elseif xprime(h)<infac_ll*if0 | xprime(h)>infac_hl*if0
        break;
      end

% Calculates cost
      cost_new = 0;                                        % Resets cost value to sum over datasets
      for i_set = 1:num_dataset
        if (strcmp(ptgpstr,'cubic') | strncmp(ptgpstr,'o',1) | strncmp(ptgpstr,'t',1))
          Vnew = {zeros(1,5);[0 0 0 0 xprime(1) 0 0 0 5*xprime(1)]; ...
                  [zeros(1,6) xprime(2) 0 0 0 -21*xprime(2) 0 0]};
        else
          for i_s = 1:size(V,2)                            % i_s indexes sites.
            vd(i_s,find(vd(i_s,:))) = xprime(1:nV(i_s)); 
            Vnew(:,1) = {vd(i_s,1:5);vd(i_s,6:14);vd(i_s,15:27)};
          end
        end
        [sptmp,pk] = saficf_genspec(J,T(i_set),Vnew,xdat{i_set},Ei(i_set),freq(i_set),lineshape{i_set});
        spec_new{i_set} = elas_peak{i_set} + sptmp'.*x(n); % Add elastic and inelastic peaks
	% Calculates the cost of the difference peak positions and of the guess CF transition energies.
        pk(find(pk(:,2)<1e-2),:) = []; 
        xpk = pk(:,1); omt = [];
        for i = 1:length(xpk); 
          omt(:,i) = xpk(i)-xpeaks{i_set}; 
        end; 
        costp = sqrt(sum(min(omt.^2)));
                   
        if costflag(i_set)                                 % Calculates new cost 
          cost_new = cost_new + sqrt( sum( (spec_new{i_set} - ydat{i_set}).^2 ) ) + costp;
        else
          cost_new = cost_new + sqrt( sum( (spec_new{i_set} - ydat{i_set}).^2 ./ (edat{i_set}.^2) ) ) + costp;
        end
      end
% Acceptance test
      if (cost_new<cost) | (rand < exp((cost-cost_new)/Tsa)) 
        x = xprime;
        cost = cost_new;
        n_u(h) = n_u(h) + 1;
      end
      if (cost<cost_opt)
        x_opt = x; 
        V_opt = Vnew
        cost_opt = cost;
      end;

    end                                                    % Step 4: add 1 to h; if h<=n goto 1
   end                                                     %         else h=0; add 1 to j

% ------------------------------------------------------------------------------------------------- %
% Step 5: Updates step vector
% ------------------------------------------------------------------------------------------------- %
   for i_u = 1:n
     if n_u(i_u) > (0.6*Ns)
       vprime(i_u) = v(i_u) * (1 + ci(i_u)*(n_u(i_u)/Ns - 0.6)/0.4 );
     elseif n_u(i_u) < (0.4*Ns)
       vprime(i_u) = v(i_u) / (1 + ci(i_u)*(0.4 - n_u(i_u)/Ns)/0.4 );
     else
       vprime(i_u) = v(i_u);
     end
   end
   
   v = vprime;
   n_u = zeros(1,n);

   %cost_opt
   %toc
  end                                                      % Step 5: set j=0; add 1 to m;

% ------------------------------------------------------------------------------------------------- %
% Step 6: Reduce temperature
% ------------------------------------------------------------------------------------------------- %
  Tsa = rT * Tsa
  cost_opt
  f_ustar(k+Nepsilon+1) = cost;

  for i_set = 1:num_dataset
    if ~isempty(handle) && ishandle(handle(i_set))  % Updates graphs
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
  
  toc
% ------------------------------------------------------------------------------------------------- %
% Step 7: Determine terminating criterion
% ------------------------------------------------------------------------------------------------- %
  term_fl = 0;
  for i_e = 1:Nepsilon
    if f_ustar(k+Nepsilon+1)-f_ustar(k+Nepsilon+1-i_e) < epsilon
      term_fl = term_fl + 1;
    end
  end
  if (f_ustar(k+Nepsilon+1)-cost_opt < epsilon) && term_fl==Nepsilon
    endex = 1
    break;                                                 % Ends search
  else
    x = x_opt;
    V = V_opt;
    cost = cost_opt;
  end

end                                                        % Step 6: set m=0; add 1 to k

% ------------------------------------------------------------------------------------------------- %
% Post-SA processing 
% ------------------------------------------------------------------------------------------------- %

% Draws final graphs
for i_set = 1:num_dataset
  if ~isempty(handle) && ishandle(handle(i_set))
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
