function mask = Mask_Vortex(P, order, orientation, x, y) 

  nx = cos(orientation*2*pi/360.0);
  ny = sin(orientation*2*pi/360.0);

  mask = abs( exp(1i*(x*nx + y*ny)/P) + exp(1i*order*atan2(x,y))  )^2;

end
