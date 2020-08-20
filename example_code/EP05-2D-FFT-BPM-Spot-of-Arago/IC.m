function ic = IC(x)

w0 = 2.0e-03;
w1 = 5.0e-03;

ic = exp(-(x^2/w1^2)^4);
if     abs( x ) < w0
ic = 0.0;
end

end
