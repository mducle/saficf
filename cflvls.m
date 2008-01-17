function peaks = cflvls(Hcf,T,flags,nt,A)
% cflvls - Plots the CF energy levels and calculates its dipole matrix elements.
%
% Syntax:  peaks = cflvls(Hcf,flag,T,nt)
%
% Inputs:  Hcf   - a (2J+1)x(2J+1) matrix representing the CF Hamiltonian in the |J,Jz> 
%                  basis. Alternative, a cell array {J B} or {A B}, where A = [L S J] 
%                  B = {B2 B4 B6} contains the angular momentum quantum numbers and
%                  crystal field parameters of the magnetic ions. 
%          T     - a scalar = the temperature in Kelvin at which to calculate the dipole 
%                  matrix element. If this is omitted, it is automatically set to zero K.
%          flags - flags to control the plotting and output formats. 
%                  flags = [1 1] - Plots CF levels and outputs in two columns (default) 
%                  flags = [0 1] - Does not plot CF levels, still outputs in two columns
%                  flags = [0 0] - Does not plot CF levels, outputs in complex numbers.
%                                  This last case is for ease of plotting: 
%                                  e.g. plot(peaks) vs. plot(peaks(:,1),peaks(:,2))
%          nt    - either a scalar = the number of levels to plot transitions from. 
%                  or a vector with the levels to plot - e.g. nt = [1 5] to plot ground 
%                  state and 5th excited state.
%                  If omitted, the function will just plot the ground state transitions.
%          A     - vector [L S J] with the angular momentum quantum numbers of the
%                  magnetic ion to be calculated. If the last input argument is a three
%                  vector, the function will calculate the transition cross-section in
%                  barns instead of just the transition dipole matrix elements.
%
% Outputs: peaks - is either a two column vector with the energy transfer and associated
%                  dipole matrix element, or a complex number vector with the real part
%                  as the energy and the imaginery part as the dipole matrix element.

% By Duc Le - Thu Aug 10 00:33:52 BST 2006 - duc.le@ucl.ac.uk

% mdl - updated 060817: improved complex output so `plot(cflvls(Hcf,T,[0 0])` works 
%     - updated 070123: fixed mag_op_j(J) error [was mag_op_J(4)] and added Hcf checks
% mdl - updated 070126: added transition cross-section (in Barns) output.

% This file is part of the SAfiCF package. Licenced under the GNU GPL v2 or later. 

% Physical constants. Taken from NIST Reference on Constants, Units, and Uncertainty, 
% http://physics.nist.gov/cuu/Constants/
mu_B  = 927.400949e-26;    % J/Tesla - Bohr magneton
k_B   = 1.3806505e-23;     % J/K     - Boltzmann constant
Q_e   = 1.60217653e-19;    % C       - Charge of electron
k_Be  = k_B/(Q_e/1000);    % meV/K   - Boltzmann constant
N_A   = 6.0221415e23;      % Avogadro's number
r0    = 0.5291772108;      % pm/100  - Borh radius.
%r0 = .53822671531422537349; % For comparison with FOCUS and cfield by Peter Fabi.

% Checks for optional arguments and flags
if ~exist('Hcf')
  error('You must supply either a CF Hamiltonian matrix or a cell array {J B}.');
elseif iscell(Hcf)
  A = Hcf{1}; B = Hcf{2};
  if isscalar(A)
    J = A; Hcf = cf_hmltn(J,B);
  elseif length(A)==3;
    L = A(1); S = A(2); J = A(3);
    gJ = 1.5 + (S*(S+1) - L*(L+1)) / (2*J*(J+1));
    Hcf = cf_hmltn(J,B);
  end
elseif isnumeric(Hcf)
% Checks that CF Hamiltonian is square and hermitian
  n = size(Hcf,1);
  if n ~= size(Hcf,2)
    error('CF Hamiltonian is not square!');
  elseif n<2
    error('CF Hamiltonian too small');
  end
  J = (n-1)/2;
  if sum(sum(Hcf-Hcf')) > 1e-10
    error('CF Hamiltonian is not hermitian!');
  end
else
  error(['The first argument must be a square matrix or cell array {A B} or {J B}' ...
         'where A is a 3-vector, B = {B2 B4 B6} are the CF parameters, J is a scalar']);
end

% Checks optional input arguments
if ~exist('flags'); flags = [0 1]; elseif isvector(flags) & length(flags)==3; Ap = flags; end
if ~exist('nt');    nt    = 1;     elseif isvector(nt)    & length(nt)==3;    Ap = nt;    end
if ~exist('T');     T     = 1e-10; elseif isvector(T)     & length(T)==3;     Ap = T;     end

if exist('A') && isvector(A) && length(A)==3
  gJ = 1.5 + (A(2)*(A(2)+1) - A(1)*(A(1)+1)) / (2*A(3)*(A(3)+1));
elseif exist('Ap')  
  gJ = 1.5 + (Ap(2)*(Ap(2)+1) - Ap(1)*(Ap(1)+1)) / (2*Ap(3)*(Ap(3)+1)); 
end

% Diagonalises the CF Hamiltonian and cleans up the eigenvectors and eigenvalues.
[V,E] = eig(Hcf);
V(find(abs(V)<3e-3)) = 0;
E = diag(E) - min(min(E));

% Calculates the magnetic operators in the |J,Jz> basis
Jmat = mag_op_j(J);
Jx = Jmat(:,:,1);   
Jy = Jmat(:,:,2); 
Jz = Jmat(:,:,3);  %               2       2          2       2 
Jp = Jmat(:,:,4);  % J+   N.B. |J |  + |J |  = 2( |J |  + |J |  )
Jm = Jmat(:,:,5);  % J-          +       -          x       y

% Calculates the dipole matrix elements (prop. to scattered inelastic neutron intensity)
Trans = ( (V'*Jx*V).*conj(V'*Jx*V) + (V'*Jy*V).*conj(V'*Jy*V) + (V'*Jz*V).*conj(V'*Jz*V) ) .* (2/3);

if exist('gJ')
  Trans = ( 4*pi*(r0*gJ/2)^2 ) .* Trans;
end

% Calculates the partition function.
Z = sum(exp(-E./(k_Be*T)));

% Calculates the transition matrix elements at finite temperatures.
for ind = 1:size(E,1)
  Trans(:,ind) = Trans(:,ind) .* (exp(-E(ind)./(k_Be*T))/Z);
end

% Outputs the position (in energy) and intensity of inelastic neutron peaks.
for ind_i = 1:size(E,1)
  for ind_j = 1:size(E,1)
    ind = (ind_i-1)*size(E,1) + ind_j;
    peaks(ind,1) = E(ind_j) - E(ind_i);
  end
end
peaks(:,2) = Trans(:);

% Eliminates repeated entries due to degenerate levels.
peaks = peaks(find(abs(peaks(:,2))>1e-10),:);
for ind_i = 1:size(peaks,1)
  for ind_j = 1:size(peaks,1)
    if(ind_i~=ind_j) 
      if (abs(peaks(ind_i,:)-peaks(ind_j,:))<1e-10)
        peaks(ind_i,:) = [0 0];
        %break;
      end
    end
  end
end 
% Eliminates degenerate level transitions with different intensities.
for ind_i = 1:size(peaks,1)
  for ind_j = 1:size(peaks,1)
    if(ind_i~=ind_j) 
      if (abs(peaks(ind_i,1)-peaks(ind_j,1))<1e-10) && (abs(peaks(ind_i,2)-peaks(ind_j,2))>1e-1) 
        peaks(ind_j,2) = peaks(ind_i,2) + peaks(ind_j,2);
        peaks(ind_i,:) = [0 0];
        %break;
      end
    end
  end
end
peaks = peaks(find(peaks(:,1)),:);

% Determine which type of output format is desired.
if ~flags(2)
  peakn = [];
  peaks = peaks(:,1) + i.*peaks(:,2);
  for ind = 1:size(peaks,1)
    peakn = [peakn;real(peaks(ind));peaks(ind);real(peaks(ind))];
  end
  peaks = peakn;
end

% Plots the CF levels.
if flags(1)
  plot([0 2]',[E E],'k');
  xlim([0 2]);
  hold all;

  % Plots the levels and label them with their wavefunctions.
  jind = -4:4;
  for ind = 1:size(E,1)
    txt = ''; 
    g = V(:,ind);
    for iind = find(g)'
      txt = [txt sprintf('%0.2g|',g(iind)) int2str(jind(iind)) '\rangle +']; 
    end 
    if ind~=1 && abs(E(ind) - E(ind-1)) < 1e-2
      text(1.2,E(ind)-0.7,txt);
    else
      text(1.2,E(ind)+0.7,txt);
    end
  end

  % Cleans up the transition matrix elements.
  Trans(find(abs(Trans)<1e-2))=0;

  % Determines which transitions to plot.
  if ~isscalar(nt)
    io = nt;
  else
    io = 1:nt;
  end

  % Plots inelastic neutron transitions along with the dipole matrix element.
  num_trans = size(find(Trans(:,io)),1);
  i_t = 0;
  for iout = io
    tgnd = Trans(:,iout); 
    ind_tgnd = find(tgnd);
    for ind = ind_tgnd'
      i_t = i_t + 1;
      plot([i_t/num_trans i_t/num_trans],[E(iout) E(ind)],'^-b')
      if ind~=1 && abs(E(ind) - E(ind-1)) < 1e-2
        text(i_t/num_trans+0.01,E(ind)-0.7,sprintf('%0.2g',tgnd(ind)));
      elseif ind==1
        text(i_t/num_trans+0.01,E(ind)-0.7,sprintf('%0.2g',tgnd(ind)));
      else
        text(i_t/num_trans+0.01,E(ind)-0.7,sprintf('%0.2g',tgnd(ind)));
      end
    end
  end

  % Cleans up plot and add labels.
  set(gca,'xtick',[])
  ylabel('Energy Levels (meV)');
  hold off;
end
