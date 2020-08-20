%% 
clear all; close all;

%% load input parameters
Parameters

%% Definition of grid
dx      = LX/(NX+1);                 
index_x = 0:NX-1;
cx      = ((index_x - floor(NX/2))*dx); 

%% Initial condition = Gaussian propagating at angle

E0 = zeros(NX,1);
for x=1:NX
	E0(x) = GaussianBeam1DRotated( cx(x)+ cx(NX)*0.5, 0.0, w0, k0, rangle);
end;

Eold = E0;

%% this will hold the map of amplitude in z-x space
map = zeros(100,NX);
map(1,:) = E0;


%% PROPAGATE
tic;
zcoor = 0.0;
for k=1:zstop/zstep

  DZ = k0*zstep;

  ac = 1/4.0 + sqrt(12.0 + DZ*DZ)/(8.0*sqrt(3.0));
  bc = 1/8.0*DZ;
  Enew = CrankNicolson(k0, dx, ac, bc, Eold, NX);  
  Eold = Enew;

  ac = 1/4.0 - sqrt(12.0 + DZ*DZ)/(8.0*sqrt(3.0));
  bc = 1/8.0*DZ;
  Enew = CrankNicolson(k0, dx, ac, bc, Eold, NX);  
  Eold = Enew*exp(1i*DZ);



  map(k,:) = Eold;
  zcoor = zcoor + zstep;
end
toc;

E1 = zeros(NX,1);
for x=1:NX
	E1(x) = GaussianBeam1DRotated( cx(x) + cx(NX)*0.5, zcoor, w0, k0, rangle);
end


figure(1);
%pcolor(log10(abs(map).^2 + 1.0e-20)); shading flat; colorbar;
pcolor(abs(map).^2); shading flat; colorbar;



figure(2);
plot(cx,real(E0),'r');
hold on;
plot(cx,real(Eold),'b');
plot(cx,real(E1),'g');

figure(3);
plot(cx,real(Eold-E1.'),'b');
