%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% Linear computation domain 
%% Using tri-diagonal solver
%% Using Hadley`s ABCs
%% Including Kerr effect
%% Demonstrating spatial soliton solutions of different orders


clear all; 
close all;

%% load input parameters
Parameters_angle

%% Initial condition
%E0 = GaussianBeam1D(cx,0,w0,k0,f,kx*k0);
E0 = SpatialSoliton(cx,0,w0,k0,0.0,kx*k0);

%% nonlinear coupling
nlcoeff = 1i*k0*n2*Iunit*dz;
display(nlcoeff);

%% see how many solitons: 
solitonnumber = sqrt(k0^2*n2*Iunit*w0^2);
display(solitonnumber);

%%
Observer_Init

%% PROPAGATE
z       =  0;
Eold    = E0;
Eoldold = E0;

for k=1:Kmax
  for m=1:M
     z = z + dz;
     Enew    = CrankNicolson(dx, k0, Eold, Eoldold, NX, dz, nlcoeff);
     Eoldold = Eold;
     Eold    = Enew;
 end
  Observer_Report;
end

%%% END PROPAGATE    
    
Etarget = SpatialSoliton(cx,z,w0,k0,f,kx*k0);

figure(3)
plot(cx, abs(Eold));
hold on;
plot(cx,real(Etarget),'r');
plot(cx,abs(Etarget),'b');
plot(cx,real(Eold),'g');

figure(4)
plot(cx,real(Etarget-Eold));
