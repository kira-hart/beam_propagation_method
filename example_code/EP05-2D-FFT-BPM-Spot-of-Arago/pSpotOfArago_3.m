% 2-D beam propagator, FFT based method
% physics = spot of Arago
% this version switches off the boundary guard: look for consequences
% attempt to eliminate artficats with PBC + "infinite" beam

% parameters 
lambda = 633.0e-09;          
LX     = 1.0e-02;   % smaller domain means better resolution
NX     = 4096;
dz     = 0.04;    % propagate to the distance equal to a
LZ     = 0.1;

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;
stps = LZ/dz;
%stps = 1;           % execute only single step to avoid "compounding"
t = cputime;

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

% calculate maximal propagation angle
maxtan = kx(NX/2)/k0;
display(maxtan);


% amplitude holder, define an initial condition
w0 = 2.0e-3;
am0 = ones(NX,NX);
for x=1:NX
  for y=1:NX
    if(cx(x)^2 + cx(y)^2 < w0^2) am0(x,y) = 0.0; end
  end
end

% propagator
pxy = zeros(NX,NX);
for x=1:NX
  for y=1:NX
  %  pxy(x,y) = exp(-1i*(kx(x)^2 + kx(y)^2)/(2*k0)*dz );
    pxy(x,y) = exp(+1i*sqrt(k0^2 - (kx(x)^2 + kx(y)^2))*dz );
  end
end

% absorbing boundary guard
bxy = zeros(NX,NX);
for x=1:NX
  for y=1:NX
    bxy(x,y) = exp(-( (cx(x)^2 + cx(y)^2)/(LX*LX/5)).^8 );
  end
end

fprintf(1,'CPU time in initialization: %g\n',cputime - t);

% define one linear step
LinearStep = @(amplitudein,propagator) ifft2(propagator.*fft2(amplitudein));

t = cputime;
tic;

% execute steps
am1 = am0;

for s=1:stps
  fprintf(1,'%d out of %d , distance = %d\n',s,stps,s*dz);
  am1 = LinearStep(am1,pxy);
 % am1 = am1.*bxy;
 figure(1);
 hold off;
 plot(cx,abs(am1(1:NX,NX/2)),'b');

end

fprintf(1,'CPU time in steps: %g\n',cputime - t);
fprintf(1,'elapsed time in steps: %g\n',toc);

% plot result
figure(1);
hold off;
plot(cx,abs(am1(1:NX,NX/2)),'b');


% some diagnostic auxiliary plots
figure(2)
clf;
subplot(4,1,1)
plot(kx,real(pxy(1:NX,1)));
title('Propagator');
subplot(4,1,2)
plot(kx,abs(fft(am0(1:NX,NX/4))));
title('Initial spatial spectrum (abs)');
subplot(4,1,3)
plot(cx,abs(am0(1:NX,NX/2)));
title('Initial condition (abs)');
subplot(4,1,4)
hold on;
plot(cx,abs(am1(1:NX,NX/2)),'b');
title('Propagated complex amplitude (abs)');
hold off;

figure(3);
imagesc(abs(am1));

% compare resulting beam profile with the independent
% result generated in cmp.c program. 
% NOTE check respective initial conditions if not a perfect match. 
figure(4);
compdata = textread('amslice.dat');
plot(compdata(:,1),compdata(:,2));
hold on;
plot(cx,abs(am1(1:NX,NX/2)),'r');
