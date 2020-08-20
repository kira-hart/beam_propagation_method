%  transverse-spatial grid definition

dx      = LX/(NX+1);             % x-grid lattice spacing
index_x = 1:NX;
cx      = ((index_x-NX/2)*dx);   % coordinates in x-dimension


%%%%%%%%%%%%%%%%%% prepare CN matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nzmax = NX*NX*5;          % maximal number of non-zero matrix entries

rows = zeros(nzmax,1);    % these arrays will encode sparse matrix od DELTA
cols = zeros(nzmax,1);
valp = zeros(nzmax,1);


% mapping space to matrix or vector index
mindex = @(i,j) j*NX + i - NX;


% enumerate all non-zero matrix elements
count = 0;           % how many non-zeros enumerated so far
for i=1:NX  
for j=1:NX

count = count + 1;
loc = mindex(i,j);    % central stencil point location
rows(count) = loc;
cols(count) = loc;
valp(count) = -4;

if(i<NX)
%left 
count = count + 1;
nnl = mindex(i+1,j);  % nearest neighbor location 
rows(count) = loc;
cols(count) = nnl;
valp(count) = 1;
end

if(i>1)
%left 
count = count + 1;
nnl = mindex(i-1,j);
rows(count) = loc;
cols(count) = nnl;
valp(count) = 1;
end

if(j>1)
%down 
count = count + 1;
nnl = mindex(i,j-1);
rows(count) = loc;
cols(count) = nnl;
valp(count) = 1;
end

if(j<NX)
%up 
count = count + 1;
nnl = mindex(i,j+1);
rows(count) = loc;
cols(count) = nnl;
valp(count) = 1;
end;

end  % j loop
end  % i loop

% discrete Laplacian operator
DELTA = sparse(rows(1:count),cols(1:count),valp(1:count),NX*NX,NX*NX,count);

% construct Lplus, Lminus

delta   = dz/(4*k0*dx^2);   
idelta  = 1i*delta;

KD = sparse(1:NX*NX,1:NX*NX,1); % Kronecked delta (identity) matrix

LP =  KD + idelta*DELTA;    
LM =  KD - idelta*DELTA;    

