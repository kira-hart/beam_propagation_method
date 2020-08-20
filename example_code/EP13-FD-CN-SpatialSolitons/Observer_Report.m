
 QQ = (Eold);
 quant_vs_z(:,k+1)=(abs(QQ(:)));
 locations_z(k+1)=z;
          
 figure(2)
 
 subplot(2,1,1)            
 plot(cx,real(QQ(:)),'b')       
 hold on;     
 plot(cx,abs(QQ(:)),'r')   
 hold off;    
 text(0.6*max(cx),0.8*max(abs(QQ)),['z=' num2str(z) 'm'])
 xlabel('x ({\mu}m) ');ylabel('I(x) (arb.u.)')
           
 subplot(2,1,2)
 pcolor(locations_z, cx, quant_vs_z); shading flat;
 axis([0 zmax -max(cx)*0.2 max(cx)*0.2])
 xlabel('z (m)'); ylabel('radius (m)');
              
 pause(tpause)          
              
  
