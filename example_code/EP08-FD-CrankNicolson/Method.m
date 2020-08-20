%  transverse-spatial grid
  
  dx      = LX/(NX+1);             % x-grid lattice spacing
  index_x = 1:NX;
  cx      = ((index_x-NX/2)*dx)';  % coordinates in x-dimension

 % a super-Gaussian boundary guard - set to no absorption for now
  bguard  = exp(-0*(cx/(0.95*LX)).^16); 


% prepare CN matrices

  delta    = dz/(4*k0*dx^2); 
  idelta   = 1i*delta;

  LP_diag = zeros(NX,1);
  LM_diag = zeros(NX,1);
  LM_diag(1:NX) = 1 + 2*idelta;
  LP_diag(1:NX) = 1 - 2*idelta;

  offdiag = ones(NX-1,1)*idelta;
    
  LM  =  sparse(1:NX-1,2:NX,  -offdiag,  NX,NX)+...
         sparse(1:NX,  1:NX,   LM_diag,  NX,NX)+...
         sparse(2:NX,  1:NX-1,-offdiag,  NX,NX);

  LP  =  sparse(1:NX-1,2:NX,  +offdiag,  NX,NX)+...
         sparse(1:NX,  1:NX,   LP_diag,  NX,NX)+...
         sparse(2:NX,  1:NX-1,+offdiag,  NX,NX);
      
