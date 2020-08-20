%% Illustraion of BPM based on implicit finite-difference scheme
%% 1D case 
%% verification of the accuracy order: w.r.t. dx (grid spacing)

clear all; close all;

%% load input parameters
Parameters1

NX = 32;
dz = 1.0e-8;


jmax = 12;
ervals = zeros(jmax,1);
dxvals = zeros(jmax,1);

for j=1:jmax

dxvals(j) = LX/NX;

%% Definition of grids and discretization matrix
Method

%% Set the distance at which solutions will be compared
z = zmax;

%% Generate start-from data
Esource = GaussianBeam1D(cx,z,w0,k0,f,0);


%% PROPAGATE SOME STEPS
stps = 10;
Eold = Esource;
for s=1:stps
Enew = LM\Eold;   % solution to  LM Enew = Eold
Eold = Enew;
end

%%% END PROPAGATE    

%% Generate what the solution should converge to
Etarget = GaussianBeam1D(cx,z+stps*dz,w0,k0,f,0);

  
%% Evaluate error
errarray = abs((Etarget - Enew));

ervals(j) = max(abs(errarray));

fprintf(1,'%d %g %g\n',NX,dxvals(j), ervals(j));

%% Improve resolution by doubling number of grid points
NX = 2*NX;



end


plot(log10(dxvals),log10(ervals),'*');

shft =4;
partialfit = polyfit(log10(dxvals(end-shft-4:end-shft)),log10(ervals(end-shft-4:end-shft)),1)
