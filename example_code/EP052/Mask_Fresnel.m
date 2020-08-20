function mask = Mask_Fresnel(F, n, lambda, x, y)

  zone = floor( (x*x + y*y)/(lambda*F) );

  if zone > n 
   mask = 0.0;
   return;
  end

 if rem(zone,2) ~= 1
  mask = 1.0;
 else
   mask = 0.0;
 end
end
