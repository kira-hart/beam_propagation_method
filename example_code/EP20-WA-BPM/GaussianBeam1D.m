
function gb = GaussianBeam1D( x, z, wx, k, f, kx) 

  gb = 1.0;

  if wx <= 0.0 
    return;
  end

  if f ~= 0.0 
    aux1 =  1.0/(wx*wx) + i*k/(2.0*f);
  else
    aux1 =  1.0/(wx*wx);
  end

  aux2 =  1.0/(1.0 + 2.0*i*z/k*aux1);

  gb =  sqrt(aux2)*exp(-(x*x)*aux1*aux2 + i*kx*x);
end
