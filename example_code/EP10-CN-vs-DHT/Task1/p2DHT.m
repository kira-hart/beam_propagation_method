% 2-D radially symmetric beam propagator, DHT based method
% Comparison with the CN BPM solution for Poisson`s spot

% parameters - should be same as in p1FFT.m
lambda = 100.0e-06;          
LR     = 1.5e-02;
NR     = 2048;
dz     = 0.0001;
LZ     = 0.01;

% derived parameters
k0 = 2*pi/lambda;
stps = LZ/dz;

% prepare Hankel transofrm
tstart = cputime;
HT = myDHT(LR,NR);
fprintf(1,'CPU time in initialization: %g\n',cputime - tstart);


% amplitude holder, define an initial condition
am = zeros(1,NR);
for x=1:NR
    am(x) = IC(sqrt(HT.cr(x)^2));
end
am = am';

% propagator'
pr = zeros(1,NR);
for x=1:NR
  pr(x) = exp(-1i*(HT.kt(x)^2)/(2*k0)*dz );
 % pr(x) = exp(i*dz*(sqrt(k0*k0 - HT.kt(x)^2) - k0) );
end
pr = pr';


% define one linear step'
LinearStep = @(amplitudein,propagator) HT.T*( propagator.*(HT.T*amplitudein) );


tstart = cputime;

% evolve amplitude
for s=1:stps
  fprintf(1,'executing %d out of %d, distance = %d\n',s,stps,s*dz);
  am = LinearStep(am,pr);
end

fprintf(1,'CPU time in propagation: %g\n',cputime - tstart);

% symmetrize radial functions for better viewwing
am2save = abs([flipud(am);am]).^2;
rc2save = [flipud(-HT.cr);HT.cr];

plot(abs([flipud(am);am]).^2);

% save for comparison with p1FFT.m result
fout = fopen('amplitude_vs_radius_DHT.dat','w');
for x=1:2*NR
    fprintf(fout,'%g %g\n',rc2save(x),am2save(x));
end


