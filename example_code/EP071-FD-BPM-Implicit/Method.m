%  transverse-spatial grid
  
  dx      = LX/(NX+1);             % x-grid lattice spacing
  index_x = 1:NX;
  cx      = ((index_x-NX/2)*dx)';  % coordinates in x-dimension

 % a super-Gaussian boundary guard - set to no absorption for now
  bguard  = exp(-0*(cx/(0.95*LX)).^16); 


% prepare the discretization matrix

  delta    = dz/(2*k0*dx^2); 
  idelta   = 1i*delta;

  LM_diag = zeros(NX,1); %these will become too large to store
  LM_diag(1:NX) = 1 + 2*idelta; %comes from the laplacian

  offdiag = ones(NX-1,1)*idelta;
    
  LM  =  sparse(1:NX-1,2:NX,  -offdiag,  NX,NX)+... %off diagonal
         sparse(1:NX,  1:NX,   LM_diag,  NX,NX)+...%diagonal
         sparse(2:NX,  1:NX-1,-offdiag,  NX,NX); %off diagonal

      
