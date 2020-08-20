%% Illustraion of BPM based on implicit discretization scheme. 
%% 1D case 

clear all; close all

%% load input parameters
Parameters0

%% Definition of grids an Crank-Nicholson related matrices
Method

%% Initial condition 
E0 = GaussianBeam1D(cx,0,w0,k0,f,0);

%%
Observer_Init

%% PROPAGATE
%% LM below is the matrix/operator that defines the implicit method
Eold = E0;

    for k=1:Kmax
        for m=1:M
            nz=(k-1)*M+m;
            z = nz*dz;
            Enew = LM\Eold;   % solution to  LM Enew = Eold
            Eold = Enew;
        end
        
        Observer_Report

    end

%%% END PROPAGATE    
  



