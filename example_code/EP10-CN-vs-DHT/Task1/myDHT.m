
function DHT = myDHT(rmax, Nr);

DHT.rmax  = rmax;
DHT.Nr    = Nr;


% read BesselJ0 zeros from the provided file
allzeros  = textread('J0zeros.dat');

% select subset of zeros to use
subzeros  = allzeros(1:Nr);

% this represents the outer boundary
lastzero  = allzeros(Nr+1);

% create symmetric auxiliary
auxmatrix = besselj(0, 1/lastzero*subzeros*(subzeros'));

% auxiliary vector that enters eqn.(2.53):
auxvector = (besselj(1,subzeros)).^(-2);

% this is (almost) Y(m,i) in (2.53) 
gsldht    = 2/lastzero*auxmatrix*diag(auxvector);

% this should be pretty close to an identity matrix
nrmaux    = gsldht*gsldht;

normalization = 1.0/sqrt(sum(diag(nrmaux))/Nr)

differencefromone = 1- normalization

% this will hold the transformation matrix
DHT.T     = normalization*gsldht;

% coordinates (radial) in real space
DHT.cr    = rmax*subzeros/lastzero;

% trnasverse wavenumbers
DHT.kt    = subzeros/rmax;

