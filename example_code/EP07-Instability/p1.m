% 1-D beam propagator, attempt for an explicit method
% illustrates: instability takes over, 

% parameters
lambda = 800.0e-09;
wx     = 40.0e-06;
LX     = 1.0e-03;
NX     = 512;
dz     = 25.0e-06;
LZ     = pi*wx*wx/lambda;
steps  = LZ/dz;

% derived parameters
k0 = 2*pi/lambda;
dk = 2*pi/LX;
dx = LX/NX;

% coordinates
cx = dx*(linspace(0,NX-1,NX)-NX/2); 

% amplitude holder, define an initial condition
% beam propagating at angle
am0 = GaussianBeam1D(cx,0.0,wx,k0,0,k0/50.0);

% ceffcient for the FD scheme
propcoeff = -1i*dz/(2.0*k0*dx*dx);
display(propcoeff);
pi*wx*wx/lambda

% define one linear step, assuming PBC
LinearStep = @(ampin,pcoeff) ampin + pcoeff*(circshift(ampin,[0,1]) -2.0*ampin + circshift(ampin,[0,-1]));


am1 = am0;
for s=1:steps
am1 = LinearStep(am1,propcoeff);
fprintf('%g\n', s*dz); 
plot(cx,real(am1));
pause(0.1);
end



