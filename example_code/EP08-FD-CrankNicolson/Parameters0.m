%% Input beam parameters
  lambda = 800e-9;            % wavelength 
  w0     = 1.0e-03;           % beam radius (waist)
  f      = 0.5;               % beam focus

%% Computational grid etc.    
  LX = 0.02;                  % x-box size (transverse dimension)
  NX = 512;                   % number of points in x-dimension

  
%% Integration control
  zmax   =  1.0;              % propagation distance 
  Nz     =  1000;             % number of propagation steps
  dz     =  zmax/Nz;          % propagation step (cm)

%% Report (diagnostics) control
  M      = 30;
  Kmax   = Nz/M + 1;

%% Derived parameters
  k0     = 2*pi/lambda;      % wavenumber in vacuum 
  Ldf    = k0*w0^2/2;        % the Rayleigh (diffraction) length  


  
 
