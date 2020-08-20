%% Illustraion of BPM based on Crank-Nicolson Scheme. 
%% 2D case 

clear all; close all

%% load input parameters
Parameters3

%% Definition of grids an Crank-Nicholson related matrices
tic;
Method
inittime = toc;
fprintf(1,'time in init: %g\n',inittime);

%% Initial condition 
E0 = zeros(NX,NX);
for i=1:NX    
for j=1:NX    
E0(i,j) = GaussianBeam2D(cx(i),cx(j),0,w0,k0,f);
end
end


%% transition to vector representation
Eold = reshape(E0,NX*NX,1);

figure(2);

%% PROPAGATE
%% LM and LP below are matrices/operators that define C-N method

z=0;
for k=1:Nz
z = z + dz;
B = LP*Eold;            % obtain rhs = LM Eold
Enew = LM\B;            % solution to  LM Enew = B
Eold = Enew;

fprintf(1,'step %d z=%g\n',k,z);
plot(real(Eold( NX*NX/2+1:NX*NX/2+NX ))); 
pause(0.01);
if z >= zstop 
break;
end
end
 
%%% END PROPAGATE    
  
Enew = reshape(Eold,NX,NX);

%% Analytic "traget" 
E1 = zeros(NX,NX);
for i=1:NX    
for j=1:NX    
E1(i,j) = GaussianBeam2D(cx(i),cx(j),z,w0,k0,f);
end
end


figure(1);
plot(cx,real(Enew(:,NX/2)));
hold on;
plot(cx,real(E1(:,NX/2)),'r');
plot(cx,real(E0(:,NX/2)),'g');

figure(2);
imagesc(real(Enew));


