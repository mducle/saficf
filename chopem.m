function res = chopem(Et,Ei,freq,chopper,Instr,phi,x2)
% chopem - calculates resolution as function of energy transfer for chopper spectrometers
%
% Syntax:  res = chopem(Et,Ei,freq,chopper,Instr)
%
% Inputs:  Et      - vector - energy transfer to calculate resolution at.
%          Ei      - scalar - incident energy
%          freq    - scalar - frequency of chopper
%          chopper - string - type of chopper: 
%                           - 's' - sloppy chopper
%                           - 'a' - a chopper
%                           - 'b' - b chopper
%                           - 'c' - c chopper
%          Instr   - string - Instrument: 'HET', 'MARI', 'MAPS', 'MERLIN'
%          phi     - scalar - the scattering angle in degrees
%          x2      - scalar - moderator to detector distance in metres
%
% Outputs: res     - vector - resolution at specified energy transfer.

% Duc Le - Mon Aug 21 21:26:40 BST 2006 - duc.le@ucl.ac.uk

% This is a port of the CHOP program by T.G. Perring, using some parts of the MCHOP port
% by J.W. Taylor for use in the SAFICF package to estimate inelastic peak widths.

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
mu_B  = 927.400949e-26;    % J/Tesla - Bohr magneton
k_B   = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e   = 1.60217653e-19;    % C - Charge of electron
N_A   = 6.0221415e23;      % Avogadro's number
hbar  = 1.05457168e-34;    % Js - planck's constant over 2*pi
m_e   = 9.1093826e-31;     % kg - mass of electron
m_n   = 1.67492728e-27;    % kg - mass of neutron

% Catches errors
if ~exist('Et') || ~isnumeric(Et)
  error('You must supply a vector of the energy transfer at which to calculate resolution.');
end
if ~exist('Ei') || ~isscalar(Ei)
  error('You must supply the incident energy as a scalar.');
end

%------------------------------  Default Parameters  ---------------------------------%

if ~exist('Instr')   Instr   = 'HET'; end
if ~exist('chopper') chopper = 's';   end
if ~exist('freq')    freq    = 200;   end   % Chopper frequency (Hz)
if ~exist('phi')     phi     = 19;    end   % Scattering angle (degrees)
if ~exist('x2')      x2      = 2.512; end   % For HET chose the 2.5m banks by default

%----------------------------  Instrument Parameters  --------------------------------%

% Sets lengths l1, l2, l3; moderator params; default V sample depending on instrument.
% Reference: LENGTHS.FOR, MODERATOR.FOR and SAMPLE.FOR from CHOP source
% NB HET and MAPS uses the water moderator, and parameters shown are correct,
%    MARI uses the methane moderator, and parameters shown are NOT correct!
if isstr(Instr)
  if strcmp(upper(Instr),'HET')
    x0 = 10.00;      % Moderator-Chopper distance (m)
    xa = 7.19;       % Aperture-Chopper distance (m)
   %x2 = 2.512;      % 2.5m bank on HET rather than 4m bank
   %x2 = 4.042;      % 4m bank on HET
    x1 = 1.82;       % Chopper-Sample distance (m)
    imod = 2;        % Water moderator - use Chi-squared 2 pulse shape
    delta0 = 38.6;   % Characteristic thickness of moderator (mm)
    deltaG = 0.5226; %   Characteristic energy dependent thickness in mm/sqrt(meV)
    wa = 66.667e-3;  % Aperture full width in (m)
    ha = 66.667e-3;  % Aperture full height (m)
    thetam = 26.7;   % Moderator angle (degrees)
    isam = 0;        % Default Vanadium sample type is flat plate
    sx = 2;          % Sample width (mm)
    sy = 40;         % Sample height (mm)
    sz = 40;         % Sample depth (mm)
  elseif strcmp(upper(Instr),'MARI')
    x0 = 10.00;      % Moderator-Chopper distance (m)                  % Incorrect
    xa = 7.19;       % Aperture-Chopper distance (m)                   % Incorrect
    x2 = 4.022;      % Sample-Detector distance (m)
    x1 = 1.739;      % Chopper-Sample distance (m)
    imod = 2;        % Water moderator - use Chi-squared 2 pulse shape % Incorrect
    delta0 = 38.6;   % Characteristic thickness of moderator (mm)      % Incorrect 
    deltaG = 0.5226; %   Char. energy-dep thickness / mm/sqrt(meV)     % Incorrect
    B1 = 0;          %   Inverse time constant for moderator E<130meV  % Incorrect
    B2 = 0;          %   Inverse time constant for moderator E>130meV  % Incorrect
    Emod = 0;        %   Swap over in energy from chi-sqr. to storage  % Incorrect
    wa = 66.667e-3;  % Aperture full width (m)                         % Incorrect
    ha = 66.667e-3;  % Aperture full height (m)                        % Incorrect
    thetam = 13.0;   % Moderator angle (degrees)                       % Incorrect
    isam = 2;        % Default Vanadium sample type is hollow cylinder
    sx = 20;         % Sample external radius (mm)
    sy = 19;         % Sample internal radius (mm)
    sz = 50;         % Sample height (mm)
  elseif strcmp(upper(Instr),'MAPS')
    x0 = 10.1;       % Moderator-Chopper distance (m)
    xa = 8.11;       % Aperture-Chopper distance (m)
    x2 = 6;          % Sample-Detector distance (m)
    x1 = 1.9;        % Chopper-Sample distance (m)
    imod = 2;        % Water moderator - use Chi-squared 2 pulse shape
    delta0 = 38.6;   % Characteristic thickness of moderator (mm)
    deltaG = 0.5226; %   Characteristic energy dependent thickness in mm/sqrt(meV)
    wa = 70.13e-3;   % Aperture full width (metres)
    ha = 70.13e-3;   % Aperture full height (metres)
    thetam = 32.0;   % Moderator angle (degrees)
    isam = 0;        % Default Vanadium sample type is flat plate
    sx = 2;          % Sample width (mm)
    sy = 50;         % Sample height (mm)
    sz = 50;         % Sample depth (mm)
  elseif strcmp(upper(Instr),'MERLIN')
    error('MERLIN: Sorry, not yet implemented.');
  else
    error('The Instr argument must be either ''HET'', ''MARI'', ''MAPS'', or ''MERLIN''');
  end
else
  error('The Instr argument must be a string.');
end

% Checks frequency is sensible.
if mod(freq,50) ~= 0
  error('Frequency must be multiple of 50Hz.');
end
if strcmp(upper(Instr),'MAPS') && freq<150
  error('Frequency must be greater than 150Hz on HET and MARI.');
end
if freq>600
  error('Frequency must be less than 600Hz.');
end

% Converts frequency to radians per second
omega = freq * 2*pi;

%------------------------------  Chopper Parameters  ---------------------------------%

% Sets chopper parameters
% Reference: CHOPPER.FOR from CHOP source
if strcmp(upper(Instr),'HET') | strcmp(upper(Instr),'MAR')
  switch lower(chopper)
    case 's'
%   Slit width(mm) | Slat th(mm)  | Chop rad.(mm) | Slit curv(mm) | Jitter time(us) 
      pslit = 2.28;  dslat = 0.55;   radius = 49;    rho = 1300;     tjit = 0;
    case 'a'
      pslit = 0.76;  dslat = 0.55;   radius = 49;    rho = 1300;     tjit = 0;
    case 'r'
      pslit = 1.14;  dslat = 0.55;   radius = 49;    rho = 2150;     tjit = 0;
    case 'd'
      pslit = 1.52;  dslat = 0.55;   radius = 49;    rho = 410;      tjit = 0;
  end
elseif strcmp(upper(Instr),'HET')
  switch lower(chopper)
    case 'b'
      pslit = 1.29;  dslat = 0.55;   radius = 49;    rho = 920;      tjit = 0;
    case 'c'
      pslit = 1.71;  dslat = 0.55;   radius = 49;    rho = 580;      tjit = 0;
  end
elseif strcmp(upper(Instr),'MARI')
  switch lower(chopper)
    case 'b'
      pslit = 1.14;  dslat = 0.55;   radius = 49;    rho = 820;      tjit = 0;
    case 'c'
      pslit = 1.52;  dslat = 0.55;   radius = 49;    rho = 580;      tjit = 0;
  end
elseif strcmp(upper(Instr),'MAPS')
  switch lower(chopper)
    case 'a'
      pslit = 1.087; dslat = 0.534;  radius = 49;    rho = 1300;     tjit = 0;
    case 'b'
      pslit = 1.812; dslat = 0.534;  radius = 49;    rho = 920;      tjit = 0;
    case 's'
      pslit = 2.899; dslat = 0.534;  radius = 49;    rho = 1300;     tjit = 0;
  end
else
  error('Chopper type not recognised. Please use ''a'',''b'',''c'',''d'',''s'',''r''');
end

%-----------------------------  Detector Parameters  ---------------------------------%

% Default detector parameters for HET, MARI and MAPS 
% Reference DETECTOR.FOR and VAN_VAR.FOR from CHOP source
idet = 1;  % He3 Position sensitive detectors
wd = 1;    % Detector width (m)
hd = 1;    % Detector height (m)
dd = 25;   % Detector thickness (mm)
tbin = 0;  % Binning time (microseconds)

%-------------------------------------------------------------------------------------%

% Converts values into SI units.
delta0 = delta0/1000;        % Characteristic thickness of moderator (m)
deltaG = deltaG/1000;        % in m/sqrt(meV)
thetam = thetam*pi/180;      % Moderator angle (radians)
sx = sx/1000;                % Sample dimension (m)
sy = sy/1000;                % Sample dimension (m)
sz = sz/1000;                % Sample dimension (m)
pslit = pslit/1000;          % Slit width (m)
dslat = pslit + dslat/1000;  % Slat thickness (m)
radius = radius/1000;        % Chopper radius (m)
rho = rho/1000;              % Slit radius of curvature (m)
tjit = tjit*1e6;             % Jitter time (s)
dd = dd/1000;                % Detector thickness (m)
tbin = tbin.*1e6;            % Binning time (s)
phi = phi.*pi/180;           % Scattering angle (radians)

%-------------------------  Calculates Vanadium widths  ------------------------------%

% Calculates neutron velocities for incident and final energies.
velEi = sqrt(2/m_n*Q_e/1000) * sqrt(Ei);
velEf = sqrt(2/m_n*Q_e/1000) * sqrt(Ei-Et);

%---------------------------  Moderator contribution  --------------------------------%

% Calculates the variance of the moderator pulse for various pulse shape
% Reference: subroutine TIKEDA, TCHI, and TCHI_2 in MOD_FUNS.FOR from CHOP source
if imod == 0        % Methane moderator?
  SIG = sqrt( (delta0^2) + ((deltaG^2*81.8048) / Ei) );
  A = 4.37392e-4 * SIG * sqrt(Ei);
  if (Ei(j) > 130.0);
    B = B2;
  else
    B = B1;
  end                                                                                                           
  R = exp(-Ei/Emod);
  Dtm2 = (3/(A^2)) + (R*(2-R))/(B^2);
  Dtm2 = Dtm2*1.0e-12;                   % Converts variance from microsecs^2 to secs^2
elseif imod == 1    % ?? Moderator 
  Dtm2 = ( (delta0/1.96)/ velEi )^2;  % Variance in secs^2
elseif imod == 2    % 300K Water Moderator 
  Dtm2 = ( ( (delta0+deltaG*sqrt(Ei))/1.96) / velEi )^2;   % Variance in secs^2
else
  error('Moderator type not recognised');
end

%----------------------------  Chopper contribution  ---------------------------------%

% Calculates the variance of the time pulse through the chopper
% Reference: subroutine TCHOP in CHOP_FUNS.FOR from CHOP source 
tch = x1/velEi;   % Time of flight of neutron from chopper to sample
% gamma is a parameter to determine the regime in which to calculate.
gamma = ( 2*(radius^2)/pslit ) * abs(1/rho - 2*omega/velEi);
if gamma >= 4
  Dtch2 = 0;
else
  if gamma <= 1
    gsqr = (1-gamma^4/10) / (1-gamma^2/6);
  else
    gsqr = 0.6*gamma * ((sqrt(gamma)-2)^2) * (sqrt(gamma)+8)/(sqrt(gamma)+4);
  end
  Dtch2=( (pslit/(2*radius*omega))^2 / 6) * gsqr;
end

%-----------------------------  Sample contribution  ---------------------------------%

% Calculates the sample variances in the x-, y- and z-directions
% Reference: function SAM0 in SAMPLE.FOR from CHOP source
if (isam == 0)           % Plate Sample
  v_x = sx^2 / 12;         % x-variance (m^2)
  v_y = sy^2 / 12;         % y-variance (m^2)
  v_z = sz^2 / 12;         % z-variance (m^2)
elseif (isam == 1)       % Ellipsoidal sample
  v_x = sx^2 / 20;         % x-variance (m^2)
  v_y = sy^2 / 20;         % y-variance (m^2)
  v_z = sz^2 / 20;         % z-variance (m^2)
elseif (isam == 2)       % Hollow Cylinder
  v_x = (sx^2 + sy^2) / 4; % x-variance (m^2)
  v_y = (sx^2 + sy^2) / 4; % y-variance (m^2)
  v_z = sz^2 / 12;         % z-variance (m^2)
elseif (isam == 3)       % Spherical sample
  v_x = sx^2 / 20;         % x-variance (m^2)
  v_y = sx^2 / 20;         % y-variance (m^2)
  v_z = sx^2 / 20;         % z-variance (m^2)
elseif (isam == 4)       % Solid cylinder
  v_x = sx^2 / 4;          % x-variance (m^2)
  v_y = sx^2 / 4;          % y-variance (m^2)
  v_z = sz^2 / 12;         % z-variance (m^2)
end

v_xy = 0;    % Diagonal variance?

%---------------------------  Detector Contribution  ---------------------------------%

% Calculates the detector variance
% Reference: function DETECT2 and DETECT_HE in DETECTOR.FOR from CHOP source
if (idet == 0)      % He3 PSD binned horizontally
  error('Binned He3 detectors not yet implemented.');
elseif (idet == 1)  % He3 PSD singly
  rad = dd / 2;        % Radius (m)
  atms = 10;           % Pressure (atmospheres)
  t2rad = 0.063;       % ???
  %don't call this function at the mo and approximate some values for it
  %[effic,delta,ddsqr,v_dd,v_d]=detect_he(wf,rad,atms,t2rad)
  effic=.5;
  delta=0;
  ddsqr=0.25;
  v_dd=0.01;
  v_d=0.01;
elseif (idet == 2)  % Davidson type scintillator
  error('Davidson type detectors not yet implemented.');
else
  error('Detector type not recognised and/or implemented.');
end
sigd  = wd / sqrt(12);  % Variance in width (m)
sigdz = hd / sqrt(12);  % Variance in height (m)
sigdd = sqrt(v_dd);     % Variance in depth (m)

%-------------------------------------------------------------------------------------%

rat = (velEi./velEf).^3; % Ratio of incident and scattered velocites
v_mod = Dtm2;            % Variance due to moderator (s^2)
v_ch = Dtch2;            % Variance due to chopper (s^2)
v_jit = tjit^2;          % Variance due to jitter (s^2)
tanthm = tan(thetam);    % Theta_m = moderator angle (radians)

am   = -(x1+rat*x2) / x0;
ach  = 1 + (x1+rat*x2) / x0;
g1   = 1 - omega*(x0+x1)*tanthm / velEi; 
g2   = 1 - omega*(x0-xa)*tanthm / velEi;
f1   = 1 + (x1/x0)*g1;
f2   = 1 + (x1/x0)*g2;
gg1  = g1 / ( omega*(xa+x1) );
gg2  = g2 / ( omega*(xa+x1) );
ff1  = f1 / ( omega*(xa+x1) );
ff2  = f2 / ( omega*(xa+x1) );
aa   = ( cos(gamma)/velEi - cos(gamma-phi)./velEf ) - ff2*sin(gamma);
bb   = (-sin(gamma)/velEi + sin(gamma-phi)./velEf ) - ff2*cos(gamma);
aya  = ff1 + (rat*x2/x0)*gg1;
ax   = aa  - (rat*x2/x0)*gg2*sin(gamma);
ay   = bb  - (rat*x2/x0)*gg2*cos(gamma);
a_dd = 1./velEf;

v_van_m  = am.^2   .* v_mod;       % Moderator contribution to Vanadium width (s^2) 
v_van_ch = ach.^2  .* v_ch;        % Chopper contribution to Vanadium width (s^2) 
v_van_jit= ach.^2  .* v_jit;       % Jitter contribution to Vanadium width (s^2) 
v_van_ya = aya.^2  .* (wa^2/12);   % Apperture width contribution to Vanadium width (s^2) 
v_van_x  = ax.^2   .* v_x;         % Sample x-direction contribution to Vanadium width (s^2) 
v_van_y  = ay.^2   .* v_y;         % Sample y-direction contribution to Vanadium width (s^2) 
v_van_xy = ax.*ay  .* v_xy;        % Sample diagonal contribution to Vanadium width (s^2) 
v_van_dd = a_dd.^2 .* v_dd;        % Detector contribution to Vanadium width (s^2) 

% Calculates total Vanadium width (s^2)
%v_van = v_van_m + v_van_ch + v_van_jit + v_van_ya + v_van_x + v_van_y  + v_van_xy  + v_van_dd;
% Because we have to calculated the detector contribution properly, including it here will
%   completely skew the results...
v_van = v_van_m + v_van_ch + v_van_jit + v_van_ya + v_van_x + v_van_y + v_van_xy;

% Calculates the resolution function
res = 8.747832e-4 .* sqrt((Ei-Et).^3) .* ( sqrt(v_van + (tbin^2/12)) .* 1e6 ) ./ x2;

% Converts from std. dev. to FWHM
res = res .* (2*sqrt(2*log(2)));

% If Ei < higher peaks get complex FWHM - need to be eliminated
res(find(imag(res))) = eps;

%-----------------------------  Resolution functions  --------------------------------%
function [effic,delta,ddsqr,v_dd,v_d] = detect_he(wf,rad,atms,t2rad)
% TODO: Implement Me!
