
function gb = GaussianBeam2D( x, y, z, wx, k, f) 

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

gb =  aux2*exp(-(x*x+y*y)*aux1*aux2);
end
