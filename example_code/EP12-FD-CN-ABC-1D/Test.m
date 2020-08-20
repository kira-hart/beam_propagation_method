%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% Linear domain 
%% Using tri-diagonal solver
%% Using Hadley`s ABC

clear all; 
close all;

%% load input parameters
Parameters

%% Initial condition
EB = GaussianBeam1D(cx,     0,w0,  k0,f,kx*k0);

EA = GaussianBeam1D(cx-LX/4,0,w0,  k0,f,kx*k0/2.0);


plot(cx,real(EA),'b');
hold on;
plot(cx,real(EB),'r');

pause(5.0);

E0 = 1.0*EA + 1.0*EB;
%E0 = EA + EB;


%%
Observer_Init

%% PROPAGATE
Eold = E0;
z    = 0.0;
for k=1:Kmax
  for m=1:M
     z = z + dz;
     Enew = CrankNicolson(dx, k0, Eold, NX, dz);
     Eold = Enew;
  end
  Observer_Report;
end

%%% END PROPAGATE    
    



