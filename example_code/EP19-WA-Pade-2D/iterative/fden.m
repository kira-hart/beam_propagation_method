


function r = fnum( XX )

Parameters1;
dx    = LX/(NX+1);   
coeff = 4.0/(k0*k0*dx*dx);
ikzdz = 1i*k0*dz/2.0;

p0 = reshape(XX,NX,NX);
p1 = coeff*del2(p0);
p2 = coeff*del2(p1);

pr = (p2*0.0625 + p1*0.75 +  p0) - ikzdz*(p2*0.25   + p1*0.50);

pr(:,1)  = 0;
pr(:,NX) = 0;
pr(1,:)  = 0;
pr(NX,:) = 0;

r  = reshape(pr, NX*NX, 1);



