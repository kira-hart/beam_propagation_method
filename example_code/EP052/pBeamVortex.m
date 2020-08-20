% 2-D beam propagator, FFT based method
% physics = vortex beams

% parameters 
lambda     = 633.0e-09;          
w0          = 1.5e-03;
maskorder   = 1;
maskperiod  = w0/40;

LX     = 1.5e-02;
NX     = 4096/2;
LZ     = 2.0;
dz     = 0.05;
stps   = LZ/dz;

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;

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

% amplitude holder, define an initial condition
am0 = zeros(NX,NX);
for x=1:NX
  for y=1:NX
	  am0(x,y) = exp( -(cx(x)^2 + cx(y)^2)/w0^2 )*Mask_Vortex(maskperiod, maskorder, 45, cx(x), cx(y));
  end
end

% propagator, paraxial should be sufficient 
pxy = zeros(NX,NX);
for x=1:NX
  for y=1:NX
	  pxy(x,y) = exp(-1i*(kx(x)^2 + kx(y)^2)/(2*k0)*dz );
  end
end

figure(1);
imagesc(real(am0)); colorbar; title('Initial condition, real part');

fprintf(1,'CPU time in initialization: %g\n',cputime - t);

% define one linear step
LinearStep = @(amplitudein,propagator) ifft2(propagator.*fft2(amplitudein));

t = cputime;
tic;

% execute one step
am1 = am0;

for s=1:stps
  fprintf(1,'%d out of %d , distance = %d\n',s,stps,s*dz);
  am1 = LinearStep(am1,pxy);

% show intermediate results
  figure(2);
  hold off;
imagesc(abs(am1)); colorbar;

end

fprintf(1,'CPU time in steps: %g\n',cputime - t);
fprintf(1,'elapsed time in steps: %g\n',toc);


figure(3)
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

figure(4);
imagesc(real(am1));

