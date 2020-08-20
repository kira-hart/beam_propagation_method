%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% 1D case 

clear all; close all

%% load input parameters
Parameters1

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
            Eold = Enew;            % no boundary guard for now
        end
        
        Observer_Report

    end

%%% END PROPAGATE    
  

Etarget = zeros(NX,1);
for i=1:NX    
Etarget(i) = GaussianBeam1D(cx(i),z,w0,k0,f,0);
end



figure(3)
plot(cx,real(Etarget),'r');
hold on;
plot(cx,real(Eold));

figure(4)
plot(cx,abs(real(Etarget)-real(Eold)));

% save results
fid = fopen('results.dat','w');
for i=1:NX   
    fprintf(fid,'%g %g %g\n', cx(i), real(Eold(i)), real(Etarget(i)));
end

fclose(fid);

