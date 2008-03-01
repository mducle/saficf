function [x, std]=saficf_salsqr(xdat,ydat,edat,pin,dpin,func,fcp,range)
% saficf_salsqr - Fits an arbitrary function using simulated annealing with (optional) constraints
%
% Syntax:  function [p, std]=saficf_salsqr(x,y,err,pin,dpin,func,fcp,rnge)
%
% Inputs:  xdat - vector - input independent variable
%          ydat - vector - input dependent variable (must be same size as x)
%          edat - vector - (optional) input error in y
%          pin  - vector - initial parameters
%          dpin - vector - initial changes in the parameters to use (must be same dimensions as pin)
%          func - handle - handle to the function to be fitted
%          fcp  - vector - [fraction_to_multiply_dpin_with max_number_iterations cost_func_tolerance]
%			   (optional) defaults to [0.1 20 0.0001]
%          rnge - vector - (optional) the range [-rnge...+rnge] for each parameter over which the
%                          fitted function is bounded. If not specified, the function attempts to 
%                          calculated a suitble range, which gives a cost function always less than
%                          maxCostF times the value of sqrt( sum(ydat.^2 ./ edat.^2) ).
%                          The value of maxCostF can be changed in the function m-file.
%
% Outputs: p    - vector - fitted parameters
%	   std  - vector - estimated errors in fitted parameters
%
% Ref: A. Corana, M. Marchesi, C. Martini, S. Ridella, ACM Trans. Math. Software, v13, p262-280, 1987

% Duc Le - Tue Aug 28 14:27:17 BST 2007 - duc.le@ucl.ac.uk
% This file is part of the SAFiCF package, licenced under the Gnu GPL v2. 

% update - Thu Jan 10 16:02:02 GMT 2008 - bugfixes, and finally added error estimation. 

% ----------------------------------  Some default parameters  ------------------------------------ %
maxTstep = 100;   % maximum number of temperature steps to take
maxCostF = 50;    % factor to multiply sqrt(sum(ydat.^2./edat.^2)) to find the maximum allowed cost.
maxRngeF = 500;   % factor to multiply initial value of parameter by to get the maximum range.
% ------------------------------------------------------------------------------------------------- %

% ------------------------------------------------------------------------------------------------- %

% Resets the random number generator
rand('state',sum(100*clock));

% Parses input parameters
if ~isvector(xdat) || ~isvector(ydat) || length(xdat)~=length(ydat)
  error('x, y must be vectors of the same lengths.');
elseif ~isvector(edat) || length(xdat)~=length(edat)
  if isvector(pin) && length(edat)==length(pin)
    if ~isa(dpin,'function_handle')
      error('func must be a handle to a function');
    elseif exist('func') 
      if ~isvector(func) || length(func)~=3
        error('fcp must be a vector of length 3');
      else
        fcp = func;
      end
    end
    func = dpin;
    dpin = pin;
    pin  = edat;
    edat = ones(size(xdat));
    if nargout>1
      warning(sprintf(['You have not specified error data, but want estimates of error in parameters.' ...
               '\n\tAs this functions uses the error data to calculate the covariance matrix, the' ...
               '\n\terror estimates may be too large if your actual error data is small or vice versa.']));
    end
  end
elseif ~isvector(pin) || ~isvector(dpin) || length(pin)~=length(dpin)
  error('pin and dpin must be vectors of the same lengths.');
elseif ~isa(func,'function_handle')
  error('func must be a handle to a function');
elseif exist('fcp') && (~isvector(fcp) || length(fcp)~=3)
  error('fcp must be a vector of length 3');
end

if ~exist('fcp')
  dp = -dpin.*0.1;    niter = 20;     stol = 0.0001;
else
  dp = -dpin.*fcp(1); niter = fcp(2); stol = fcp(3);
end

% Converts from syntax of speclsqr to that of Corana et al.
x = pin(:);
v = dpin(:);
n = length(pin);

% So that vectorial +/- agree
xdat = xdat(:);
ydat = ydat(:);
edat = edat(:);

% Find out which parameters we need to vary and which we don't!
variPars = find(v~=0);

% Some Simulated Annealing Parameters, after Corana et. al. 
Ns = 20;             % Number of passes to ensure step sizes give acceptance rate of approx 50%
NT = min([20,5*n]);  % Number of steps at each temperature to ensure thermal equilibrium.
Nepsilon = 4;        % Number of successive temperature reductions before testing for stopping fit.
ci = ones(1,n).*2;   % Step varying criterion intial value
rT = 0.85;           % Temperature reduction coefficient
epsilon = 1e-4;      % Convergence criterion

% Determine the range of the parameters
if ~exist('range','var')
  maxcost = maxCostF * sqrt( sum(ydat.^2 ./ edat.^2) );
  range = zeros(size(x));
  for ipar = variPars(:)'
    ht = ones(size(x)); ht(ipar) = maxRngeF*ht(ipar);
    xprime = x .* ht;   cost = maxcost*2;
    while(cost>maxcost)
      f = feval(func,xdat,xprime);
      cost = sqrt( sum( (f - ydat).^2 ./ (edat.^2) ) );
      xprime(ipar) = xprime(ipar)/1.5;
    end
    range(ipar) = xprime(ipar);
  end
else
  range = range(:);
end

% Defines the start SA temperature to allow at least 50% acceptance
% Generates 10 random steps and finds cost of each. Then takes mean of cost, and set
%   initial temperature so that at least half of steps are accepted.
cost_opt =sqrt(sum((feval(func,xdat,x)-ydat).^2./edat.^2));% Optimized parameters and cost function so far.
for ind_conf = 1:10
  xprime = x + ( ((rand(size(v))-0.5).*2).*range.*v );
  f = feval(func,xdat,xprime);
  cost(ind_conf) = sqrt( sum( (f - ydat).^2 ./ (edat.^2) ) );
  if cost(ind_conf)<cost_opt
    cost_opt = cost(ind_conf); x_opt = xprime;             % Use this cost/parameters if it's the lowest so far.
  end
end
mean_cost = sum(cost)/10;
par_steps = sum(range)/length(range)/NT;                   % Average change in parameters / step.
cost_factor = par_steps * 25 / mean_cost;
[T0,Tsa] = deal( mean_cost * cost_factor / 25 );           % 1/log((1-n)/n) = 25 for n = 0.49

% More Simulated Annealing Parameters, after Corana et. al. 
n_u = zeros(n,1);                                          % u=1,...,n; n_u is a step variation count.
f_ustar = ones(Nepsilon,1) .* cost_opt;                    % u=0,-1,...,-Nepsilon+1; stores cost of each T step.

% ------------------------------------------------------------------------------------------------- %
% Steps 1-4                                                                                         %
% ------------------------------------------------------------------------------------------------- %
for k = 0:maxTstep
  %tic
  for m = 1:NT                                             % One m-loop takes ~?? for 5 gaussians...
   for j = 1:Ns
    for h = variPars(:)'
      eh = zeros(n,1); eh(h) = 1;
      xprime = x + ( (rand-.5)*2*v(h) .* eh );             % Step 1
      if h<n
        if xprime(h) < -range(h) | xprime(h) > range(h)
          break;                                           % Step 2
        end
      elseif xprime(h)<-abs(range(h)) | xprime(h)>abs(range(h))
        break;
      end
      % Calculates cost
      f = feval(func,xdat,xprime);
      cost_new = sqrt( sum( (f - ydat).^2 ./ (edat.^2) ) );

      % Acceptance test
      if (cost_new<cost) | (rand < exp((cost-cost_new)/Tsa)) 
        x = xprime;
        cost = cost_new;
        n_u(h) = n_u(h) + 1;
      end
      if (cost<cost_opt)
        x_opt = x; 
        cost_opt = cost;
      end;

    end  % for h = variPars                                % Step 4: add 1 to h; if h<=n goto 1
   end  % for j = 1:Ns                                     %         else h=0; add 1 to j

% ------------------------------------------------------------------------------------------------- %
% Step 5: Updates step vector
% ------------------------------------------------------------------------------------------------- %
   for i_u = variPars(:)'
     if n_u(i_u) > (0.6*Ns)
       vprime(i_u) = v(i_u) * (1 + ci(i_u)*(n_u(i_u)/Ns - 0.6)/0.4 );
     elseif n_u(i_u) < (0.4*Ns)
       vprime(i_u) = v(i_u) / (1 + ci(i_u)*(0.4 - n_u(i_u)/Ns)/0.4 );
     else
       vprime(i_u) = v(i_u);
     end
   end
   
   v = vprime(:);
   n_u = zeros(1,n);

   %cost_opt
   %toc
  end  % for m = 1:Nt                                      % Step 5: set j=0; add 1 to m;

% ------------------------------------------------------------------------------------------------- %
% Step 6: Reduce temperature
% ------------------------------------------------------------------------------------------------- %
  Tsa = rT * Tsa;
  cost_opt;
  f_ustar(k+Nepsilon+1) = cost;

  %toc
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
    break;                                                 % Ends search
  else
    x = x_opt;
    cost = cost_opt;
  end

end  % for k = 1:maxTstep                                  % Step 6: set m=0; add 1 to k

% ------------------------------------------------------------------------------------------------- %
% Post-SA processing 
% ------------------------------------------------------------------------------------------------- %

f=feval(func,xdat,x_opt);
m=length(ydat);
wt=1./edat(:);
% Error estimation - lifted from speclsqr.m     % Ref: Y. Bard, Non-Linear Parameter Estimation
jac=feval('specdfdp',xdat,f,x,dp,func);
jac=jac(:,find(dp));                            % use only fitted parameters
Qinv=diag(wt.*wt);
Q=diag((0*wt+1)./(wt.^2));
resid=ydat-f;                                   %un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 
covr=diag(covr);                                %for compact storage
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);
stdresid=resid./sqrt(diag(Vy));

jtgjinv=pinv(jac'*Qinv*jac);
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
for k=1:n,
  for j=k:n,
    corp(k,j)=covp(k,j)/sqrt(abs(covp(k,k)*covp(j,j)));
    corp(j,k)=corp(k,j);
  end;
end;

std=sqrt(diag(covp));

j=1;
sig=zeros(size(x));
for i=1:length(std)
	while dp(j)==0
		j=j+1;
	end
	sig(j)=std(i);
	j=j+1;
end
std=sig;

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*jtgjinv



