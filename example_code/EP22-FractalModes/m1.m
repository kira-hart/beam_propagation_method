clear all
% 2-D beam propagator, FFT based method



% parameters 
lambda = 1.0;     % fix wavelength
aa     = 205.0*lambda;
aa     = 205.0*lambda*sqrt(5.0);
aa     = 165.0*lambda;

R1     = 1000.0*lambda;
R2     = 1500.0*lambda;
M      = R2/R1;
LZ     = (R2 - R1)/2.0;

LX     = 1.1*M*(2*aa);

NX     = 2048*1;   % assume computer can`t handle more points

%w0     = 1*10^-5;  % circular screen radius
%w1     = 5*w0;     % beam size
%dz     = 0.5*w0;   % step = one half of radius

% derived parameters
k0   = 2*pi/lambda;
dk   = 2*pi/LX;
dx   = LX/NX;


FN   = (M-1)^2*aa^2/(2*lambda*LZ);

display(M);
display(FN);


% coordinates
cx = dx*(linspace(0,NX-1,NX)-NX/2); 



al  = 2.0/3.0*pi;
Tr  = [[cos(al) sin(al)];[-sin(al) cos(al)]];
r1  = [ 0 ; 1];
r2 = Tr*r1;
r3 = Tr*r2;

% screen holder, define an initial condition

smallm = zeros(NX,NX);
for x=1:NX
  for y=1:NX
    r = [cx(x);cx(y)];
    if(  (r'*r1 > -0.5*aa)&&(r'*r2 > -0.5*aa)&&(r'*r3 > -0.5*aa)  ) %'
      smallm(x,y) = exp(+1i*k0*(cx(x)^2+cx(y)^2)/(2*R1));
    end;
  end
end

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

dz = LZ;
%px =exp(-1i*dz*(kx.^2)/(2*k0) );
%px =exp(+1i*dz*(sqrt(k0^2-kx.^2) - k0) );

% precalculate propagator
paraxial  = false;
dz = LZ;
pxy = zeros(NX,NX);
if paraxial == true;
  for x=1:NX
    for y=1:NX
      pxy(x,y) =exp(-1i*dz*(kx(x)^2+kx(y)^2)/(2*k0) );
    end
  end
else
  for x=1:NX
    for y=1:NX
	angle = sqrt(kx(x)^2 + kx(y)^2)/k0;
	angle0 = aa/R1;
        pxy(x,y) =exp(+1i*dz*(sqrt(k0^2-kx(x)^2-kx(y)^2) - k0) );
        if( angle > angle0 ) 
	  pxy(x,y) = pxy(x,y)*exp(-((angle-angle0)/(0.5*angle0))^2);
       end;
    end
  end
end

% define one linear step
LinearStep = @(amplitudein,propagator) ifft2(propagator.*fft2(amplitudein));

am0 = smallm;

%am0(NX/2+20,NX/2+54) = 1.0;

old = am0;
for i=1:100
am0 = LinearStep(am0,pxy);

am0 = am0.*largem;
am0 = LinearStep(am0,pxy);

imagesc(abs(am0).^2); colormap gray;
pause(0.01);

am0 = am0.*smallm;

aux = max(max(abs(am0)));
am0 = am0/aux;

err = max(max(abs(am0) - abs(old)));

display(i);
display(err);

old = am0;

end;

old = abs(am0).^2;
dlmwrite('matrix.dat',old,' ');



