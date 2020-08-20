% 1-D beam propagator, FFT based method
% Attempt to propagate into far field:
% This script utilizes a smaller compuational domain


% parameters
lambda = 800.0e-09;

w0     = 3.0e-05;  % slit width/2

% let us attempt a smaller lattice
LX     = 10.0e-03;
NX     = 1024*4;

% choose step and how far to propagate, show Fresnel number
dz      = 50.0e-05;
zstop   = 0.04;
fresnel = w0^2/(zstop*lambda)

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;

% coordinates
cx = dx*(linspace(0,NX-1,NX)-NX/2); 

% transverse wavenumbers
kx = zeros(1,NX);
for k=0:NX/2
	kx(1+k) = dk*k;
end
for k=NX/2+1:NX-1
	kx(1+k) = dk*(k - NX);
end

% paraxial propagator
px = exp(-1i*(kx.*kx)/(2*k0)*dz );

%%%%%%%%%%%%%%%%%%%%%%% boundary guard %%%%%%%%%%%%%%%%%%%%%

% boundary guard: experiment wihe different shapes

bg    = ones(1,NX);

edge  = 0.25;       % portion od the domain on each side used as a guard
shape = 2.0;        % sharper transition with higher value of shape
strength = 500.0;   % something like absorption coefficient of the guard
strength = 0;       % guard not used in ths code

% right edge of the domain
x1 = (1.0-edge)*NX;
for x=x1:NX
  y = (x-x1)/(NX-x1);
  bg(x)=cos(y*pi/2.0)^shape;
end;

% left edge of the domain
x1 = edge*NX;
for x=1:x1
  y = (x1-x)/x1;
  bg(x)=cos(y*pi/2.0)^shape;
end;

% incorporate step so that the guard action does not depend on it
bg = exp(-dz*(1.0-bg)*strength);

%%%%%%%%%%%%%%%%%%% end boundary guard shapes %%%%%%%%%%%%%


% amplitude holder, define an initial condition = illuminated slit
% a smooth edge is meant to suppress high wave numbers
% experiment with the sharpnes and zoom-in on results...

am0 = exp( -((cx - cx(NX/2))/w0).^128 );

% define one linear FFT-BPM step
LinearStep = @(amplitudein,propagator) ifft(propagator.*fft(amplitudein));

distance = 0.0;
am1      = am0;

for s=0:1000
am1 = LinearStep(am1,px);
am1 = am1.*bg;             % apply boundary guard after each step
distance = distance + dz;
if(distance >= zstop); break; end;
end

% plot result
figure(1)
clf;
subplot(4,1,1)
plot(kx,real(px));
title('Propagator');
subplot(4,1,2)
plot(kx,abs(fft(am0)));
title('Initial spatial spectrum (abs)');
subplot(4,1,3)
plot(cx,abs(am0));
title('Initial condition (abs)');
subplot(4,1,4)
hold on;
plot(cx,abs(am1),'r');
%plot(cx,real(am1),'b');
%plot(cx,imag(am1),'g');
title('Propagated complex amplitude (abs)');
hold off;


