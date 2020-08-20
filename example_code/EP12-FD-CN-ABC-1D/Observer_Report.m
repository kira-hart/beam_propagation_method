
 QQ = (Eold);
 mapholder(:,k+1)=log(abs(QQ(:))+1e-16);
 locations_z(k+1)=z;
          
 figure(2)
 
 subplot(2,1,1)            
 plot(cx,real(QQ(:)),'b')            
 text(0.6*max(cx),0.8*max(abs(QQ)),['z=' num2str(z) 'm'])
 xlabel('transverse location ({\mu}m) ',  'FontSize' ,12);
 ylabel('I(x) (arb.u.)','FontSize' ,12)
           
 subplot(2,1,2)
 pcolor(locations_z, cx, mapholder); shading flat;
 axis([0 zmax -max(cx) max(cx)*1.0])
 xlabel('propagation distance (m)','FontSize' ,12); 
 ylabel('transverse location (m)', 'FontSize' ,12);
              
 pause(tpause)          
              
  
