
function gb = SpatialSoliton( x, z, wx, k, f, kx) 

  beta = 1/(2*k*wx^2) - kx^2/(2*k);

  gb = sech((x  - z*kx/k)/wx).*exp(1i*kx*x)*exp(1i*beta*z);

end
