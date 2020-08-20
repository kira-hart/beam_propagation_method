function ic = IC(x)

w0 = 2.0e-03;
w1 = 8.0e-03;
fr = w0/20.0;

fr = 0.0;

if (x<w0-fr)
  ic = 0.0;
end;

if ((w0-fr<x)&&(x<w0))
  ic = ( sin( pi/2.0*(x - w0 + fr)/fr  ) )^2.0;
end;

if (x>w0)
  ic = exp(-(x^2/w1^2)^8.0);
end;

end
