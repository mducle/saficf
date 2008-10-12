function Hz = zeeman_hmltn_dir(A,H,Hdir)
% zeeman_hmltn_dir - calculates the directional Zeeman Hamiltonian 
% given by: H = g u  J.H
%            z     B - -  
% where H is the applied magnetic field, u_B is the Bohr magneton
% and g is the lande g-factor.
%
% Syntax:  Hz = zeeman_hmltn(A,H,Hdir) 
%
% Inputs:  A = [L S J] with L,S,J being the angular momentum quantum 
%              numbers of the ground state multiplet.
%          H is the magnitude of the magnetic field in tesla, in the 
%              form H = [H_min:step:H_max]
%          Hdir = [Hx Hy Hz] is a vector describing the direction H.
%
% Outputs: Hz is a (2J+1)x(2J+1) matrix representing the Zeeman Hamiltonian

% By Duc Le (2005) - duc.le@ucl.ac.uk

% Physical constants. Taken from G. Woan, The Cambridge Handbook of 
% Physics Formulas, CUP 2000
mu_B = 5.78838263e-2;    % meV/T - Bohr magneton

% Makes equations look nicer:
L = A(1); S = A(2); J = A(3);

%                                    3   S(S+1) - L(L+1)
% Calculates the Lande g-factor: g = - + ---------------
%                                    2       2J(J+1)
g = 1.5 + (S*(S+1) - L*(L+1)) / (2*J*(J+1));

% Normalises Hdir
Hdir = Hdir ./ sqrt(Hdir * Hdir');

% Calculates the matrix elements <Vi|J|Vj>
Jmat = mag_op_j(J);
Jmat = Jmat(:,:,1).*Hdir(1) + Jmat(:,:,2).*Hdir(2) + Jmat(:,:,3).*Hdir(3);

% Calculates the Hamiltonian.  Given by: -g u  J.H
%                                            B - - 
index = 0;
for h = H
  index = index + 1;
  for ind_i = 1:(2*J+1)
    for ind_j = 1:(2*J+1)
      Hz(ind_i,ind_j,index) = -g * mu_B * Jmat(ind_i,ind_j) * h;
    end
  end
end 
