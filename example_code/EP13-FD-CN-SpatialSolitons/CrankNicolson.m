function Enew = CrankNicolson(dr, k0, Eold, Eoldold, n, zstep, nlcoef)  

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
  A(row) =      - idelta*A(row);
  B(row) = +1.0 - idelta*B(row);
  C(row) =      - idelta*C(row);
end


% Bare implementation of ABCs: without taking precautions
% This is the place to modify when ABCs cause problems!

% Absorbing boundary conditions - LHS modification
B(n) = B(n) + C(n)*Eold(n)/Eold(n-1);
B(1) = B(1) + A(1)*Eold(1)/Eold(2);

% Absorbing boundary conditions - ghost values to be used in RHS
LBNDRY = 0.0;
RBNDRY = 0.0;
LBNDRY = Eold(1)/Eold(2)*Eold(1);
RBNDRY = Eold(n)/Eold(n-1)*Eold(n);

%% Define the right-hand-side of the linear systsem to solve, using boundary ghost values:

row = 1;
   R(row) = Eold(row) + idelta*( LBNDRY      -2.0*Eold(row) + Eold(row+1)  );

for row=2:n-1
   R(row) = Eold(row) + idelta*( Eold(row-1) -2.0*Eold(row) + Eold(row+1)  );
end

row = n;
   R(row) = Eold(row) + idelta*( Eold(row-1) -2.0*Eold(row) + RBNDRY  );

%% RHS modifcation due to nonlinearity
for row=1:n
  R(row) = R(row) +  nlcoef*3/2*Eold(row)*abs(Eold(row))^2 -  nlcoef*1/2*Eoldold(row)*abs(Eoldold(row))^2;
end


%% Solve the linear system L_minus E_new = L_plus E_old

Enew = TDMAsolver(A, B, C, R);

