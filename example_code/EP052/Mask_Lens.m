function mask = Mask_Lens(F, k0, x, y) 
  mask = exp(-1i*(x*x + y*y)*k0/(2*F));
end



