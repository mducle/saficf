function out = threej(A)
% threej - Calculates the value of the 3-j symbol
%
% Syntax:  out = threej(A)
%
% Inputs:  A   - A 2x3 matrix: ( j1 j2  J )
%                              ( m1 m2 -M )
% Outputs: out - is the value of the 3-j symbol, a scalar.
%
% Based on the program THREEJ.BAS from the Crystal Field Handbook by Newman and Ng (CUP 2000)
% The 3j symbols are given in mathworld.wolfram.com

% By Duc Le - duc.le@ucl.ac.uk 2005
% mdl - updated 060810 for saficf: cleaned up some constructions. Made equations look nicer.

% This file is part of the SAfiCF package. 
% Licenced under the GNU GPL v2 or later. 

[szrow, szcol] = size(A);
if (szrow ~= 2) | (szcol ~= 3)
  error('Input requires a 2x3 matrix')
end

% To make the following equations look nicer
j1 = A(1,1); j2 = A(1,2); J = A(1,3);
m1 = A(2,1); m2 = A(2,2); M = A(2,3);

% Selection rules:  for the symbols:  ( j1 j2  J )
%                                     ( m1 m2 -M )
%
% The parameters j1,j2,J,m1,m2,-M are integers or half integers
for i = 1:6
  if A(i) ~= 0 && mod(2*A(i),1) 
    out = 0; return;
  end
end 
% m1 + m2 = M  i.e  m1 + m2 - M = 0
if sum(A(2,:)) ~= 0
  out = 0; return;
% The triangular inequality |j1-j2| <= J <= j1+j2
elseif (J < abs(j1-j2)) | (J > j1+j2)
  out = 0; return;
% m1 = -j1, ... +j1
elseif (m1 > j1) | (m1 < -j1)
  out = 0; return;
% m2 = -j2, ... +j2
elseif (m2 > j2) | (m2 < -j2)
  out = 0; return;
% M = -J, ... +J
elseif (-M > J) | (-M < -J)
  out = 0; return;
% Integer perimeter rule: j1 + j2 + J is integer
elseif mod(sum(A(1,:)),1)
  out = 0; return;
end

% The analytical expression for the Wigner 3j symbol is given by the Racah formula
% ( j1 j2 J )        (j1-j2-M)  
% ( m1 m2 M ) = (-1)^          sqrt(T(j1,j2,J)) sqrt((j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(J+M)!(J-M)!)
%                                               t
%               * SUM _____________________(-1)^________________________________
%                  t  t!(J-j2+t+m1)!(J-j1+t-m2)!(j1+j2-J-t)!(j1-t-m1)!(j2-t+m2)!
%
% where T(a,b,c) is the triangle coefficient: (a+b-c)!(a-b+c)!(-a+b+c)! /
%                                                    (a+b+c+1)!


% Calculates sum over t. The t cover all integers for which the factorials exist.
sum_t = 0;
tmin = max( 0, max( j2-m1-J, j1+m2-J ) );
tmax = min( j1+j2-J, min( j1-m1, j2+m2 ) );
for t = tmin:tmax
  sum_t = sum_t + (-1)^t / prod(factorial([t J-j2+t+m1 J-j1+t-m2 j1+j2-J-t j1-t-m1 j2-t+m2]));
end

% Calculates the triangle function
T = prod(factorial([j1+j2-J j1-j2+J -j1+j2+J])) / factorial(j1+j2+J+1);

% Calculates the 3j symbol
out = (-1)^(j1-j2-M) * sqrt(T) * sqrt(prod(factorial([j1+m1 j1-m1 j2+m2 j2-m2 J+M J-M]))) * sum_t;

