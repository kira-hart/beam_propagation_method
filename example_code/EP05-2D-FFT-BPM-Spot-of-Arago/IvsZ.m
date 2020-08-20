% 2-D beam propagator, FFT based method
% physics = spot of Arago

% parameters, physical 
lambda =   ;
w0 =  ;  %disk-size
w1 =  ;  %beam-size

% parameters, computational
LX     =   % domain size
NX     =   % grid-points
dz     =   % step
LZ     =   % total distance

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;
stps = floor(LZ/dz);

t = cputime;

% coordinates
cx = dx*(linspace(0,NX-1,NX)-NX/2+1); 

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
    rad = sqrt(cx(x)^2 + cx(y)^2);
    am0(x,y) = exp(- (rad/w1)^8.0);
    if(rad < w0)
      am0(x,y)  = 0.0;
    end
  end
end

% propagator
pxy = zeros(NX,NX);
for x=1:NX
  for y=1:NX
%    pxy(x,y) = exp(-1i*(kx(x)^2 + kx(y)^2)/(2*k0)*dz );
     pxy(x,y) = exp(+1i*dz*(sqrt(k0*k0 - kx(x)^2 - kx(y)^2)-k0));
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

centralIntensity= zeros(1,stps);
targetIntensity = zeros(1,stps);
locs            = zeros(1,stps);

for s=1:stps
  fprintf(1,'%d out of %d , distance = %d\n',s,stps,s*dz);
  am1 = LinearStep(am1,pxy);
  am1 = am1.*bxy;


centralIntensity(s) = abs(am1(NX/2,NX/2))^2;
targetIntensity(s) = (s*dz)^2/((s*dz)^2 + w0^2);
locs(s) = s*dz;

figure(5);
plot(locs(1:s),centralIntensity(1:s));
hold on;
plot(locs(1:s),targetIntensity(1:s),'or');
hold off;
xlabel('propagation distance');
ylabel('on-axis intensity');
pause(0.1);
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
%plot(compdata(:,1),compdata(:,2));
%hold on;
plot(cx,abs(am1(1:NX,NX/2)),'r');
