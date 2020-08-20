%% Input beam parameters
  lambda  = 1.06e-06; % wavelength 
  w0 = 2.0/sqrt(2.0)*2.828e-06;      % beam radius (waist)
  f  =+0.0e+00;      % beam focus
  angle = +35;       % angle in deg
  
%% Computational grid etc.    
  LX = 100e-06;             % x-box size (transverse dimension)
  NX = 8000;                % number of points in x-dimension
  
%% Derived parameters
  rangle = (angle/360.0)*2*pi;
  k0 = 2*pi/lambda;
  kx = 2*pi/lambda*sin(rangle);

%% Integration
zstop = 50.0e-06;
zstep = 0.5e-06;  

zstep = 3.0/k0;
