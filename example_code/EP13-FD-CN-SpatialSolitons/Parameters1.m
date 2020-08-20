%% Input beam parameters
  lambda = 800e-9;            % wavelength 
  w0     = 1.0e-04;           % beam radius (waist)
  f      = 0.0e+00;           % beam focus
  kx     = 0.0;               % transverse wavenumber ( relative to k0)
  Iunit  = 1.0*1.62113e+14;   % sets electric field amplitude units 
  n2     = 1.0e-20;           % nonlinear index typical of condensed media

%% Computational grid etc.    
  LX = 0.004;             % x-box size (transverse dimension)
  NX = 200;               % number of points in x-dimension

  
%% Integration control
  zmax   =  1.1234;              % propagation distance 
  Nz     =  2000;             % number of propagation steps
  dz     =  zmax/Nz;          % propagation step (cm)

%% Report (diagnostics) control
  M      = 15;
  Kmax   = Nz/M + 1;

%% Derived parameters
  k0     = 2*pi/lambda;      % wavenumber in vacuum 
  Ldf    = k0*w0^2/2;        % the Rayleigh (diffraction) length  

%% Grid
  dx      = LX/(NX+1);             % x-grid lattice spacing
  index_x = 1:NX;
  cx      = ((index_x -NX/2)*dx);  % transverse coordinates 
 
