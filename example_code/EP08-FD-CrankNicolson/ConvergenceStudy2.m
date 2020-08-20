%% Illustraion of BPM based on C-N finite-difference scheme
%% 1D case 
%% this is a convergence study in which we improve the resolution AND step
%% note that this program runs a bit longer...

clear all; close all;

%% load input parameters
Parameters1

NX = 4096*2;
NX = 64
dz = 1.0e-6;


jmax = 12;
ervals = zeros(jmax,1);
dzvals = zeros(jmax,1);

for j=1:jmax

dzvals(j) = dz;

%% Definition of grids and discretization matrix
Method

%% Set the distance at which solutions will be compared
z = zmax;

%% Generate start-from data
Esource = zeros(NX,1);
for i=1:NX    
Esource(i) = GaussianBeam1D(cx(i),z,w0,k0,f,0);
end

figure(1);
plot(cx,real(Esource));
hold on;

%% PROPAGATE SOME STEPS
stps = 2.0e-6/dz;
Eold = Esource;
for s=1:stps
    RHS  = LP*Eold; % prepare RHS
    Enew = LM\RHS;  % solution to  LM Enew = RHS
    Eold = Enew;
end

figure(1);
plot(cx,real(Eold),'r');
hold off;

%%% END PROPAGATE    

%% Generate what the solution should converge to
Etarget = zeros(NX,1);
for i=1:NX    
Etarget(i) = GaussianBeam1D(cx(i),z+stps*dz,w0,k0,f,0);
end
  
%% Evaluate error
errarray = abs((Etarget - Enew));
ervals(j) = max(abs(errarray));

%% show progress
fprintf(1,'%d %g %g\n',NX,dzvals(j), ervals(j));

%% Improve resolution by doubling number of grid points
NX = NX*2;

%% ... and decrease the integration step to one half
dz = dz/2.0;

end

figure(2);
plot(log10(dzvals),log10(ervals),'*');

%%  change the index selection to see how the slope changes:
partialfit = polyfit(log10(dzvals(8:12)),log10(ervals(8:12)),1)
