function out = lrnd(c,size_x,size_y)
% Generates a random number with a lorentzian distribution, of width 'c'.
%
% Syntax:  out = lrnd(c,size_x,size_y)
%
% Inputs:  c      - scalar - the width, or control parameter of the distribution 
%          size_x - scalar - number of columns in matrix of random numbers
%          size_y - scalar - number of rows in matrix of random numbers
%
% Outputs: out    - matrix - size_x x size_y matrix of lorentzian distributed random numbers.
%
% This function uses the inverse transform method.
%
% Reference: Luc Devroye, Non-Uniform Random Variate Generation, Chapter 2, Springer-Verlag
%   (1986). http://cg.scs.carleton.ca/~luc/rnbookindex.html

if ~exist('size_y') 
  if ~exist('size_x') 
    size_x = 1;
    size_y = 1;
  else
    size_y = size_x;
  end
end

out = rand(size_x,size_y);

% The lorentzian distribution is:  c / (c^2 + x^2);
% Applying the inverse function to uniformly distributed random numbers give us lorentzian
%   distributed random numbers by the inverse transform method
% Reference: p29, ibid.
out = c .* tan(pi.*(out-0.5))/100;
