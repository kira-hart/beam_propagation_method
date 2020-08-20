function mask = Mask_Axicon(angle, k0, x, y) 
  mask = exp(-1i*sqrt(x*x + y*y)*k0*angle);
end
