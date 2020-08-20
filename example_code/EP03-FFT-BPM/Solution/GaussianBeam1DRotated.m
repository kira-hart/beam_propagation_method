function gb = GaussianBeam1DRotated( x, z, w, k, angle)

zt = +cos(angle)*z + sin(angle)*x;
xt = -sin(angle)*z + cos(angle)*x;


aux = 1.0 + 2.0*1i*zt/(k*w^2);

gb = exp( -xt.*xt./(w^2*aux) + 1i*k*zt)./sqrt(aux);

end
