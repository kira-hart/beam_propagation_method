%% Input beam parameters
  lambda = 800e-9;            % wavelength 
  w0     = 1.0e-04;           % beam radius (waist)
  f      = 0.0e+00;           % beam focus
  kx     = -2.5e-3;           % relative to k0
  kx     = +3.0e-3;
  kx     = +10.0e-3;
%  kx     = 0.0;
  Iunit  = 1.0*1.62113e+14;
  n2     = 1.0e-20;

%% Computational grid etc.    
  LX = 0.003;             % x-box size (transverse dimension)
  NX = 250;               % number of points in x-dimension

  
%% Integration control
  zmax   =  0.3;              % propagation distance 
  Nz     =  300;              % number of propagation steps
  dz     =  zmax/Nz;          % propagation step (cm)

%% Report (diagnostics) control
  M      = 1;
  Kmax   = Nz/M + 1;

%% Derived parameters
  k0     = 2*pi/lambda;      % wavenumber in vacuum 
  Ldf    = k0*w0^2/2;        % the Rayleigh (diffraction) length  

%% Grid
  dx      = LX/(NX+1);               % x-grid lattice spacing
  index_x = 1:NX;
  cx      = ((index_x -NX/2)*dx);    % coordinates 
 
