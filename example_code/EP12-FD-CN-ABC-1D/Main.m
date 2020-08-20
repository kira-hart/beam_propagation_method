%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% Linear domain 
%% Using tri-diagonal solver
%% Using Hadley`s ABC

%% Different aprameter sets and input conditions induce
%% regime in which the ABCs work well or fail...


clear all; 
close all;

%% load input parameters. 
%% Diff/inspect the two files, comment/uncomment...
%% either:
%Parameters

%% or:
Parameters

%% Initial condition: 
%% The last parameter controls the propagation angle - this
%% is something to explore: residual reflectivity does depend on angle
E0 = GaussianBeam1D(cx,0,w0,k0,f,kx*k0);

%%
Observer_Init

%% PROPAGATE
Eold = E0;
z    = 0.0;
for k=1:Kmax
  for m=1:M
     z = z + dz;

%% This is "bare" implementation of ABC
 %    Enew = CrankNicolson(dx, k0, Eold, NX, dz);

%% This implementation is a bit more robust
 Enew = CrankNicolson(dx, k0, Eold, NX, dz);

     Eold = Enew;
  end
  Observer_Report;
end

%%% END PROPAGATE    
    



