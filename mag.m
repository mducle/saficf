function M = mag(A,H,Jdir,T,Hcf)
% outputs the magnetisation due to the crystal field parallel to a specified direction 
% 
% Syntax:  M = mag(A,H,Jdir,T,Hcf) 
%
% Inputs:  A = [L S J] with L,S,J being the angular momentum quantum 
%              numbers of the ground state multiplet; 
%          H is the applied field in teslas of the form of a vector 
%              H = [H_min:step:H_max]; 
%          Jdir is the direction of the magnetisation in form [x y z];
%          T is the temperature in Kelvins of the form of a vector 
%              T = [T_min:step:T_max];
%          Hcf is the crystal field hamiltonian which may be obtained from the
%              function xtalfld_hmltn_stev.
%
% This function calculates the magnetisation from the equation:
%                             ---      ^        -Ei(H)   / ---     -Ei(H)
% M(H,T) = g u  << J >> = gu  >   <i|J.H|i> exp(------) /  >   exp(------)
%             B             B ---i                kT   /   ---i      kT
%
% where u_B is the Bohr magneton, g is the lande g-factor, |i> are 
% eigenstates of the total Hamiltonian, H = Hcf + Hz, where Hz = -g u_B J.H
% and the hat (^) over the H indicates a unit vector. 
% Ei is the energy of the ith level, k is Boltzmann constant and T is the 
% temperature.

% By Duc Le (2006) - duc.le@ucl.ac.uk

% Derivation of the equation (from Jensen and Mackintosh, Rare Earth Magnetism)
% Starting with the Hamiltonian:
%                                H = H   - g u  J.H
%                                     cf      B
% we get the energy: Ei = <i|H  |i> - g u  <i|J.H|i> = E     - g u  <i|J.H|i>
%                             cf         B              cf,i      B     
% 
% The magnetisation is given by eqn (1.2.22):
%
%     N ---   dEi     -Ei    / ---     -Ei
% M = - >   - --- exp(---)  /  >   exp(---) = equation given above without N/V
%     V ---i   dH      kT  /   ---i     kT    factor.

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
mu_B  = 927.400949e-26;    % J/Tesla - Bohr magneton
k_B   = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e   = 1.60217653e-19;    % C - Charge of electron
N_A   = 6.0221415e23;      % Avogadro's number

% Makes equations look nicer:
L = A(1); S = A(2); J = A(3);

%                                    3   S(S+1) - L(L+1)
% Calculates the Lande g-factor: g = - + ---------------
%                                    2       2J(J+1)
g = 1.5 + (S*(S+1) - L*(L+1)) / (2*J*(J+1));

% Calculates the zeeman hamiltonian
Hz = zeeman_hmltn_dir(A,H,Jdir);

% Normalises the direction vector
Jdir = Jdir ./ sqrt(Jdir*Jdir');

% Calculates the magnetic operator matrix
Jmat = mag_op_j(J);
Jmat = Jmat(:,:,1).*Jdir(1) + Jmat(:,:,2).*Jdir(2) + Jmat(:,:,3).*Jdir(3);

% Defining beta = 1/kT here saves computation later on.
beta = 1 ./ (k_B*T);

for ind_H = 1:size(H,2)
% Calculates the total Hamiltonian as a function of field (last index)
  Hmltn = Hcf + Hz(:,:,ind_H);

% Calculates the eigenvectors V and eigenvalues (enegies) E
% Where:            ---
%        | V  >  =  >    a  |j, j    >
%           i       ---i  i      z,i
%
  [V, Edummy] = eig(Hmltn);

% Reduce the energy levels to a vector.
  Edummy = sort( Edummy(logical(eye(size(Edummy,1)))) );
% Sets energy levels relative to lowest level.
  E = Edummy - min(Edummy);

% Converts energy levels from meV to J
  E = E .* (Q_e/1000);

  mj = -J:J;
  for ind_j = 1:size(V,2)
% Calculates the matrix elements <Vi|J.H|Vi>
    me = V(:,ind_j)' * Jmat * V(:,ind_j);

% Calculates the elements of the expectation of J.H: <i|J.H|i> exp(-Ei(H)/kT)
    JH(:,ind_j) = me * exp(-beta .* E(ind_j));

% Calculates the elements of the partition function: exp(-Ei(H)/kT)
    Z(:,ind_j) = exp(-beta .* E(ind_j));
  end

  for ind_T = 1:size(T,2)
% Calculates the expectation <<J.H>> = sum(<i|J.H|i>exp(-Ei/kT)) / sum(exp(-Ei/kT))
    exp_JH(ind_H,ind_T) = sum(JH(ind_T,:)) / sum(Z(ind_T,:));
  end

end

% Calculates the magnetisation M(ind_H,ind_T) per unit volume per atom;
M = g * exp_JH;  % in u_B/atom
%M = g * mu_B * exp_JH .* N_A; % J/T/m^3/mol = A/m/mol
%To get magnetisation in emu/g (cgs units): Get magnetisation in J/T/m^3/kg by
%  M = g*exp_JH * mu_B * N_A/(molar mass). Then multiply by 4pi*10^{-7}.
%To get susceptibility in emu/mol (cgs units): Take M in J/T/m^3/mol, divide by
%  field in Tesla, and divide result by 4pi
%NB. H = B/mu_0 = B/(4pi*10^{-7}) where B is in Tesla. 
%  and chi = M/H not M/B => chi = M/(field in Tesla)*(4pi*10^{-7}) in m^3/mol.
%NB. 1T = 1Netwon/A/m 
