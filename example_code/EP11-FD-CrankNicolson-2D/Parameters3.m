%% Input beam parameters
  lambda = 800e-9;            % wavelength 
  w0     = 0.5e-03;           % beam radius (waist)
  f      = 0.4;               % beam focus

%% Computational grid etc.    
  LX = 0.005;                  % x-box size (transverse dimension)
  NX = 300;                   % number of points in x-dimension
  
%% Integration control
  zstop  =  0.4;            % propagation distance 
  Nz     =  500;            % MAX number of propagation steps
  dz     =  0.0075;          % propagation step (cm)

%% Derived parameters
  k0     = 2*pi/lambda;      % wavenumber in vacuum 
  Ldf    = k0*w0^2/2;        % the Rayleigh (diffraction) length  


  
 
