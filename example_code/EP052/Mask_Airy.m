function mask = Mask_Airy(P, x, y) 
  mask = exp(1i*2*pi*( (x/P)^3 + (y/P)^3 ) );
end
