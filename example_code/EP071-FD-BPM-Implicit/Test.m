%% Illustraion of BPM based on Finite-Difference Implicit 
%% 1D case, paraxial 

clear all; close all

%% load input parameters
ParametersCMP

%% Definition of grids and implicit evolution operator matrix
Method

%% Initial condition 
E0 = GaussianBeam1D(cx,0,w0,k0,f,0);

%%
Observer_Init

%% PROPAGATE
Eold = E0;

    for k=1:Kmax
        for m=1:M
            nz=(k-1)*M+m;
            z = nz*dz;
            Enew = LM\Eold;    % solution to  LM Enew = Eold
            Eold = Enew;
        end
        
        Observer_Report

    end

%%% END PROPAGATE    
  
Etarget = GaussianBeam1D(cx,z,w0,k0,f,0);


figure(3)
plot(cx,real(Etarget),'r');
hold on;
plot(cx,real(Eold));

figure(4)
plot(cx, abs(real(Etarget)-real(Eold)));

max(abs(Etarget-Eold))

