function Enew = CrankNicolson(k0, dx, acoeff, bcoeff, Eold, n)  

%% Define matrix L_minus tri-diagonals

A = ones(n,1);
B = -2.0*ones(n,1);
C = ones(n,1);
R = zeros(n,1);

Lcoeff = (acoeff - 1i*bcoeff)/(k0*dx*k0*dx);
Rcoeff = (acoeff + 1i*bcoeff)/(k0*dx*k0*dx);

for row=1:n
A(row) =      +Lcoeff*A(row);
B(row) = +1.0 +Lcoeff*B(row);
C(row) =      +Lcoeff*C(row);
end


%%BC
LBNDRY = 0;
RBNDRY = 0;

if abs(Eold(n-1)) > 0
%B(n) = B(n) + C(n)*Eold(n)/Eold(n-1);
%RBNDRY = Eold(n)/Eold(n-1)*Eold(n);
end;

if abs(Eold(2)) > 0
%B(1) = B(1) + A(1)*Eold(1)/Eold(2);
%LBNDRY = Eold(1)/Eold(2)*Eold(1);
end;




%% Define the right-hand-side of the linear systsem to solve, using boundary ghost values:

row = 1;
R(row) = Eold(row) + Rcoeff*( LBNDRY      -2.0*Eold(row) + Eold(row+1)  );


for row=2:n-1
R(row) = Eold(row) + Rcoeff*( Eold(row-1) -2.0*Eold(row) + Eold(row+1)  );
end

row = n; 
R(row) = Eold(row) + Rcoeff*( Eold(row-1) -2.0*Eold(row) + RBNDRY       );


%% Solve the linear system L_minus E_new = L_plus E_old

Enew = TDMAsolver(A, B, C, R);

