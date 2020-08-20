clear all
% 2-D beam propagator, FFT based method

% parameters 
lambda = 1.0;           % fix wavelength

aa     = 144.95*lambda;  % size of smaller mirror
aa     = 204.93*lambda*sqrt(2.0);  % size of smaller mirror


aa     = 504.93*lambda;  % size of smaller mirror


LX     = 3.5*aa;          % computational domain size
NX     = 1024*2;        

M      = 1.5;
LZ     = 2*1000*lambda;

R1     = 2*LZ/(M-1);
R2     = 2*M*LZ/(M-1);

% derived parameters
k0   = 2*pi/lambda;
dk   = 2*pi/LX;
dx   = LX/NX;

% coordinates
cx = dx*(linspace(0,NX-1,NX)-NX/2); 

FN   = (R2/R1-1)*aa^2/(2*lambda*LZ);
display(FN);


% set up triangular mirror
al = 0.4*pi; 
Tr = [[cos(al) sin(al)];[-sin(al) cos(al)]];

r1 = [ 0 ; 1];
r1 = Tr*r1;

al = 2.0/3.0*pi;
Tr = [[cos(al) sin(al)];[-sin(al) cos(al)]];

r2 = Tr*r1;
r3 = Tr*r2;

% small mirror aperture and phase-screen
smallm = zeros(NX,NX);
for x=1:NX
  for y=1:NX
    r = [cx(x);cx(y)];
    if(  (r'*r1 > -0.5*aa)&&(r'*r2 > -0.5*aa)&&(r'*r3 > -0.5*aa)  ) %'
      smallm(x,y) = exp(+1i*k0*(cx(x)^2+cx(y)^2)/(2*R1));
    end;
  end
end


% large mirror aperture and phase-screen
largem = zeros(NX,NX);
for x=1:NX
  for y=1:NX
    if(  cx(x)^2 + cx(y)^2 < 2.5*aa*aa  )
      largem(x,y) = exp(-1i*k0*(cx(x)^2+cx(y)^2)/(2*R2));
    end;
  end
end


% transverse wavenumbers
kx = zeros(1,NX);
for k=0:NX/2
	kx(1+k) = dk*k;
end
for k=NX/2+1:NX-1
	kx(1+k) = dk*(k - NX);
end


% precalculate propagator
for x=1:NX
  for y=1:NX
     % non-paraxial propagator
     pxy(x,y) =exp(+1i*LZ*(sqrt(k0^2-kx(x)^2-kx(y)^2) - k0) );

     % kill extreme wavenumbers
     angle  = sqrt(kx(x)^2 + kx(y)^2)/k0;
     angle0 = aa/R1;
     if( angle > angle0 ) 
      pxy(x,y) = pxy(x,y)*exp(-((angle-angle0)/(0.5*angle0))^2);
     end;

   end
end


% define one linear step
LinearStep = @(amplitudein,propagator) ifft2(propagator.*fft2(amplitudein));

% initial beam guess

am0 = rand(NX,NX);
am0 = am0.*smallm;

% auxiliary to check convergence
old = am0;

% cavity round-trip loop
for i=1:200

am0 = LinearStep(am0,pxy);
am0 = am0.*largem;
am0 = LinearStep(am0,pxy);


am0 = am0.*smallm;

aux = max(max(abs(am0)));
am0 = am0/aux;
err = max(max(abs(am0) - abs(old)));

display(i);
display(err);

old = am0;

imagesc(abs(am0).^2); %colormap gray;
pause(0.01);


end;

old = abs(am0).^2;
dlmwrite('matrix.dat',old,' ');



