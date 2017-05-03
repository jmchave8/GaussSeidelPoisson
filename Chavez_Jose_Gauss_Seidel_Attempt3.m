%Matlab Code to solve Poisson's equation with Gauss Seidel Method with the following conditions in the problem statement. 
% Jose Chavez  1161146 
clear all; clc; 
tic 
%% Given Conditions 
ax = -pi; 
ay = -pi; 
bx = pi; 
by = pi; 
N=input('Value of X Intenal Nodes='); % Number of points on the internal nodes for N and M%
M=input('Value of Y Internal Nodes='); 
Me=M+2; %Number of points including exterior boundary points for Ne and Me%
Ne=N+2; 
% this generates the x and y values that will be used to calculate 
x = linspace(-pi,pi,Ne); 
y = linspace(-pi,pi,Me); 
%% 

U = ones(Ne,Me); %U initial guess %

% For loop solving for right hand side with F equation with i,j indices% 
for i=1:length(x); 
    for j=1:length(y); 
F(i,j) = cos ( (0.5*pi)* (2*((x(i)-ax) / (bx - ax))+1 )).*sin( pi*((y(j)-ay) / (by -ay))); 
    end 
end 

%% Boundary Conditions for "top" and "bottom" 

% Bottom boundary values 
phi = ((x - ax).^2 ) .* sin( (pi *(x- ax)) / (2*(bx-ax)) ) ; 

% Top boundary values 
psy = cos (pi*(x-ax)).*cosh(bx-x); 

% place these known values in the solution grid 
U(:,1) = phi; 
U(:,Me) = psy; 
%% Left and Right Boundary points 
%   Using the given neumann condition yields special cases of the Gauss-siedel iteration that can be used along entire "side" boundaries. 
%   F and U is computed in solution grid 
% Multipliers that are used in the iterations. 
DX = 2*pi/(N+1); 
B = 1/DX.^2;
DY = 2*pi/(M+1); 
C = 1/DY.^2;
DEN = -2*(B+C); 

% Normalize Multipliers%
B = B/DEN; 
C = C/DEN; 
F = F/DEN; 
DEN = 1; 
error=10; 
error_iterations=0;
% check for diagonal dominance of elements 
abs(DEN) >= abs(2*B+2*C)
while error>10^-10; 
    W=U; 
for P = 1:1000; 

for j = 2:M+1; 
     
    % Left boundary 
    U(1,j) = DEN*(  F(1,j) - (2*B)*U(2,j) - C*U(1,j-1) - C*U(1,j+1) ); 
    % Right Boundary 
    U(N,j) = DEN*(  F(N,j) - (2*B)*U(N-1,j) - C*U(N,j-1) - C*U(N,j+1) ); 
end 

%% Gauss-Siedel iterating the general U equation%

for i = 2:N+1; 
    for j = 2:M+1; 
        U(i,j) = DEN*(  F(i,j) - C*U(i+1,j) - C*U(i-1,j)- B*U(i,j+1) - B*U(i,j-1) ); 
    end 
end 
end 
error=abs(max(max(((W-U)./W)))); 
error_iterations=error_iterations+1;
end 
toc 
error_iterations
figure 
subplot(1,2,1),surf(U),xlabel('y axis'),ylabel('x axis');
subplot(1,2,2),contour(U),xlabel('y axis'),ylabel('x axis');
