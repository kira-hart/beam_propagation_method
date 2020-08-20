function Enew = CrankNicolson(dr, k0, Eold, n, zstep)  

% This function solves:
% L_minus E_new = L_plus E_old
% where L are C-N matrices
% arguments: delta_r, k0 = 2 Pi/lambda, initial condition array, NR, zstep



%% Define matrix L_minus tri-diagonals

A = ones(n,1);
B = -2.0*ones(n,1);
C = ones(n,1);
R = zeros(n,1);

idelta = 1i*zstep/(4.0*k0*dr*dr);

for row=1:n
  A(row) = -idelta*A(row);
  B(row) = +1.0 - idelta*B(row);
  C(row) = -idelta*C(row);
end

LBNDRY = 0.0;
RBNDRY = 0.0;


if( abs(Eold(n-1)) > 1.0e-16 )

  kx = -1i*log(Eold(n)/Eold(n-1))
  if(real(kx)>0)
    B(n) = B(n) + C(n)*Eold(n)/Eold(n-1);
    RBNDRY = Eold(n)/Eold(n-1)*Eold(n);
  end;
end;

if( abs(Eold(2)) > 1.0e-16 )

  kx = -1i*log(Eold(1)/Eold(2))

  if(real(kx)>0)
    B(1) = B(1) + A(1)*Eold(1)/Eold(2);
    LBNDRY = Eold(1)/Eold(2)*Eold(1);
  end;

end;


%% Proceed to define the right-hand-side of the linear systsem to solve

row = 1;
   R(row) = Eold(row) + idelta*( LBNDRY      -2.0*Eold(row) + Eold(row+1)  );

for row=2:n-1
   R(row) = Eold(row) + idelta*( Eold(row-1) -2.0*Eold(row) + Eold(row+1)  );
end

row = n;
   R(row) = Eold(row) + idelta*( Eold(row-1) -2.0*Eold(row) + RBNDRY  );


%% Solve the linear system L_minus E_new = L_plus E_old

Enew = TDMAsolver(A, B, C, R);

