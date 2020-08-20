%% Illustraion of BPM based on implicit discretization scheme. 
%% 1D case, extracting numerical dispersion, using spatial-spatial spectrum technique

clear all; close all

%% input parameters
  lambda = 1.0e-06;            % wavelength 

%% Computational grid etc.    
  LX = 0.0004;                % x-box size (transverse dimension)
  NX = 1024;                   % number of points in x-dimension

%% Integration control
  dz     =  2.0e-07;          % propagation step 

%% Derived parameters
  k0     = 2*pi/lambda;      % wavenumber in vacuum 

%% Definition of grids an Crank-Nicholson related matrices
Method

%% Initial condition 
E0 = rand(NX,1) - 0.5 + 1i*( rand(NX,1)-0.5  );

%% holder for the spatial-spatial spectrum
NS = 1024*4;
holder = zeros(NS,NX);

%% longitudinal and transverse wavenumber coordinates
ind = 1:NS;
kz = pi/dz/NS*(ind - NS/2);
ind = 1:NX;
kx = pi/dx/NX*(ind - NX/2);

%% PROPAGATE
%% LM below is the matrix/operator that defines the implicit method
Eold = E0;

    for k=1:NS
      RHS  = LP*Eold;
      Enew = LM\RHS;   % solution to  LM Enew = RHS
      Eold = Enew;

% store snapshot of the beam profile
holder(k,:) = Eold;

% sanity check during the simulation...
figure(1);
plot(cx,real(Eold));
ylim([-0.6 0.6]);
xlabel('transverse coordinate');
hold off;

    end

%%% END PROPAGATE    
  

spectrum = fft2(holder);

figure(2);
imagesc(kz,kx,log10(abs(fftshift(spectrum))));
colorbar;
xlabel('transverse wavenumber');
ylabel('propagation constant');
ylim([-6e5 1e5]);
