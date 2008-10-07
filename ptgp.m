function allowed = ptgp(ptgp_string)
% ptgp - gives the non-zero CF parameters for a given point group.
%
% Syntax:  allowed = ptgp(ptgp_string)
%
% Inputs:  ptgp_string - a string with the point group name. 
%
% Outputs: allowed     - a cell array {B2 B4 B6} of ones and zeros of the allowed
%                        crystal field parameters.

% Reference - cfield.c file from the CFIELD program by Peter Fabi (nee Hoffmann)
%   included in the mcphase distribution, http://www.mcphase.de
%
% By Duc Le - Thu Aug 10 00:33:27 BST 2006 - duc.le@ucl.ac.uk

% This file is part of the SAfiCF package. 
% Licenced under the GNU GPL v2 or later. 

if ~exist('ptgp_string') || ~isstr(ptgp_string)
  error('You must supply the point group name as a string.');
end

% Common symmetries:
triclinic    = { [1 0 1 1 1]                    % B2 | a != b != c
                 [1 1 1 1 1 1 1 1 1]            % B4 | alpha != beta != gamma
                 [1 1 1 1 1 1 1 1 1 1 1 1 1] }; % B6 | 
monoclinic   = { [0 0 1 0 1]                    % B2 | a != b != c
                 [1 0 1 0 1 0 1 0 1]            % B4 | alpha != beta = gamma = 90
                 [1 0 1 0 1 0 1 0 1 0 1 0 1] }; % B6 |
orthorhombic = { [0 0 1 0 1]                    % B2 | a != b != c
                 [0 0 0 0 1 0 1 0 1]            % B4 | alpha = beta = gamma = 90
                 [0 0 0 0 0 0 1 0 1 0 1 0 1] }; % B6 |
tetragonalA  = { [0 0 1 0 0]                    % B2 | a = b != c
                 [0 0 0 0 1 0 0 0 1]            % B4 | alpha = beta = gamma = 90
                 [0 0 1 0 0 0 1 0 0 0 1 0 0] }; % B6 |
tetragonalB  = { [0 0 1 0 0]                    % B2 | a = b != c
                 [0 0 0 0 1 0 0 0 1]            % B4 | alpha = beta = gamma = 90
                 [0 0 0 0 0 0 1 0 0 0 1 0 0] }; % B6 |
trigonalA    = { [0 0 1 0 0]                    % B2 | a = b = c
                 [0 0 0 0 1 0 0 1 0]            % B4 | alpha != beta != gamma 
                 [1 0 0 1 0 0 1 0 0 1 0 0 1] }; % B6 |
trigonalB    = { [0 0 1 0 0]                    % B2 | a = b = c
                 [0 0 0 0 1 0 0 1 0]            % B4 | alpha != beta != gamma
                 [0 0 0 0 0 0 1 0 0 1 0 0 1] }; % B6 |
hexagonal    = { [0 0 1 0 0]                    % B2 | a = b != c
                 [0 0 0 0 1 0 0 0 0]            % B4 | alpha = beta = 120 != gamma
                 [0 0 0 0 0 0 1 0 0 0 0 0 1] }; % B6 |
cubic        = { [0 0 0 0 0]                    % B2 | a = b = c
                 [0 0 0 0 1 0 0 0 5]            % B4 | alpha = beta = gamma = 90
                 [0 0 0 0 0 0 1 0 0 0 -21 0 0] };%B6 |

% Finds the non-zero parameters for each point group.
switch lower(ptgp_string)
  case 'ci';  allowed = triclinic;     % Ci
  case 'cs';  allowed = monoclinic;    % Cs
  case 'c1';  allowed = triclinic;     % C1
  case 'c2';  allowed = monoclinic;    % C2
  case 'c3';  allowed = trigonalA;     % C3
  case 'c4';  allowed = tetragonalA;   % C4
  case 'c6';  allowed = hexagonal;     % C6
  case 'c2h'; allowed = monoclinic;    % C2h
  case 'c3h'; allowed = hexagonal;     % C3h
  case 'c4h'; allowed = tetragonalA;   % C4h
  case 'c6h'; allowed = hexagonal;     % C6h
  case 'c2v'; allowed = orthorhombic;  % C2v
  case 'c3v'; allowed = trigonalB;     % C3v
  case 'c4v'; allowed = tetragonalB;   % C4v - B20,B40,B44,B60,B64
  case 'c6v'; allowed = hexagonal;     % C6v
  case 'd2';  allowed = orthorhombic;  % D2
  case 'd3';  allowed = trigonalB;     % D3
  case 'd4';  allowed = tetragonalB;   % D4
  case 'd6';  allowed = hexagonal;     % D6
  case 'd2h'; allowed = orthorhombic;  % D2h
  case 'd3h'; allowed = hexagonal;     % D3h
  case 'd4h'; allowed = tetragonalB;   % D4h
  case 'd6h'; allowed = hexagonal;     % D6h
  case 'd2d'; allowed = tetragonalB;   % D2d
  case 'd3d'; allowed = trigonalB;     % D3d
  case 's4';  allowed = tetragonalA;   % S4
  case 's6';  allowed = trigonalA;     % S6
  case 't';   allowed = cubic;         % T   - B40,B60; B44 = 5B40; B64 = -21B60;
  case 'th';  allowed = cubic;         % Th
  case 'td';  allowed = cubic;         % Td
  case 'o';   allowed = cubic;         % O
  case 'oh';  allowed = cubic;         % Oh
% Non point group descriptions
  case 'triclinic';   allowed = triclinic;  
  case 'monoclinic';  allowed = monoclinic;
  case 'orthorhombic';allowed = orthorhombic;
  case 'hexagonal';   allowed = hexagonal;
  case 'cubic';       allowed = cubic;
% Special space groups
  case 'dhcp';   allowed = [trigonalB hexagonal];
  case 'p6mmc';  allowed = [trigonalB hexagonal];
  case 'd46h';   allowed = [trigonalB hexagonal];
  otherwise
    error('Point Group not recognised');    
end
