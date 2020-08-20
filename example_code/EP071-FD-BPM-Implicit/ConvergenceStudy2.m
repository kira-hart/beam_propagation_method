%% Illustraion of BPM based on implicit finite-difference scheme
%% 1D case 
%% verification of the accuracy order: w.r.t. dz (integration step)

clear all; close all;

%% load input parameters
Parameters1

NX = 4096*2;
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
Esource = GaussianBeam1D(cx,z,w0,k0,f,0);


figure(1);
plot(cx,real(Esource));
hold on;

%% PROPAGATE SOME STEPS
stps = 2.0e-6/dz;
Eold = Esource;
for s=1:stps
Enew = LM\Eold;   % solution to  LM Enew = Eold
Eold = Enew;
end

figure(1);
plot(cx,real(Eold),'r');
hold off;

%%% END PROPAGATE    

%% Generate what the solution should converge to
Etarget = GaussianBeam1D(cx,z+stps*dz,w0,k0,f,0);

  
%% Evaluate error
errarray = abs((Etarget - Enew));

ervals(j) = max(abs(errarray));

fprintf(1,'%d %g %g\n',NX,dzvals(j), ervals(j));

%% Improve resolution by doubling number of grid points
dz = dz/2.0;



end

figure(2);
plot(log10(dzvals),log10(ervals),'*');

shft =1;
partialfit = polyfit(log10(dzvals(1:4)),log10(ervals(1:4)),1)
