% initial beam 

E0 = zeros(NX,1);
for i=1:NX    
E0(i) = GaussianBeam1D(cx(i),0,w0,k0,f,0);
end




  