% 1-D beam propagator, FFT based method
% test vs exact Gaussian beam solution
% this case: beam propagates in the paraxial regime, but at steep angle
% to do: switch between paraxial and exact propagators

clear all;
close all;
 
% beam parameters
lambda        = 800.0e-09;
beamwaist     = 15.0e-06;
angle         = 0.35;

% computational domain parameters
LX     = 1.0e-03;
NX     = 1024*8;

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;
zR = pi*beamwaist^2/lambda;

% define real-space coordinates 
cx = dx*(linspace(0,NX-1,NX)-NX/2); 

% define transverse wavenumbers = spectral coordinates
kx = zeros(1,NX);
for k=0:NX/2
	kx(1+k) = dk*k;
end
for k=NX/2+1:NX-1
	kx(1+k) = dk*(k - NX);
end

% "integration" step
nsteps   = 20;
distance = zR;
dz       = 2*distance/nsteps;

% thi time include absolute phase in propagation

% propagator: paraxial version
px = exp(-1i*(kx.*kx)/(2*k0)*dz + 1i*dz*k0);

% nonparaxial propagator
%px = exp(+1i*dz*( (sqrt(k0*k0 - kx.*kx) - 0.0*k0) ) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial condition
am0 = GaussianBeam1DRotated(cx,-distance,beamwaist,k0,angle);

% propagated amplitude 
am1 = am0;

% define one linear step
LinearStep = @(amplitudein,propagator) ifft(propagator.*fft(amplitudein));

% execute few steps
figure(1)

currentz = 0.0;
for s=1:nsteps
am1 = LinearStep(am1,px);

currentz = currentz + dz;

plot(cx,abs(am1),'r');
hold on;
plot(cx,real(am1),'b');
hold off
pause(0.1);

end

% testing target amplitude holder, 
amt  = GaussianBeam1DRotated(cx,-distance+currentz,beamwaist,k0,angle);

figure(2)
hold on;
plot(cx,abs(am0),'b');
plot(cx,abs(amt),'r');
plot(cx,abs(am1),'b');
title('solutions compared: simulated vs exact ');
xlabel('transverse location');
ylabel('electric field envelope');
hold off;

figure(3)
hold on;
plot(cx,abs(amt-0.0*am1), ':','Color','r');
plot(cx,real(am1), 'b');
title('solutions compared:');
hold off;



