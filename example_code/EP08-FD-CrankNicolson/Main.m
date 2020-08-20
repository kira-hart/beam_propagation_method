%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% 1D case 

clear all; close all

%% load input parameters
Parameters0

%% Definition of grids an Crank-Nicholson related matrices
Method

%% Initial condition 
InitialCondition

%%
Observer_Init

%% PROPAGATE
%% LM and LP below are matrices/operators that define C-N method
Eold = E0;

    for k=1:Kmax
        for m=1:M
            nz=(k-1)*M+m;
            z = nz*dz;
            B = LP*Eold;            % obtain rhs = LM Eold
            Enew = LM\B;            % solution to  LM Enew = B
            Eold = Enew.*bguard;    % apply superGaussian apodizer
        end
        
        Observer_Report

    end

%%% END PROPAGATE    
  



