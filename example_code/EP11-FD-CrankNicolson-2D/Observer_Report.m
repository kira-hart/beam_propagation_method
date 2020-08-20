% Another example how objects (in a file) depend on what is
% defined outside - not good for a larger project


 QQ = real(Eold);
 quant_vs_z(:,k+1)=QQ(:);
 locations_z(k+1)=z;
          
 figure(2)
 
 subplot(2,1,1)            
 plot(cx,QQ(:),'b')            
 text(0.6*max(cx),0.8*max(QQ),['z=' num2str(z) 'm'])
 xlabel('x ({\mu}m) ');ylabel('I(x) (arb.u.)')
           
 subplot(2,1,2)
 pcolor(locations_z, cx, quant_vs_z); shading flat
 axis([0 zmax -max(cx)*0.3 max(cx)*0.3])
 xlabel('z (m)'); ylabel('x (m)');
              
 pause(tpause)          
              
  
