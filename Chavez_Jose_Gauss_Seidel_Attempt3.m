%Matlab Code to solve Poisson's equation with Gauss Seidel Method with the following conditions in the problem statement. 
% Jose Chavez  1161146 
clear all; clc; 

%% Given Conditions 
ax = -pi; 
ay = -pi; 
bx = pi; 
by = pi; 
N=input('Value of X Intenal Nodes='); % Number of points on the internal nodes for N and M%
M=input('Value of Y Internal Nodes='); 
tic
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
U(1,:) = phi; 
U(N+2,:) = psy; 
%% Left and Right Boundary points 
%   Using the given neumann condition yields special cases of the Gauss-siedel iteration that can be used along entire "side" boundaries. 
%   F and U is computed in solution grid 
% Multipliers that are used in the iterations. 
dx = 2*pi/(N+1); 
B = 1/dx.^2;
dy = 2*pi/(M+1); 
C = 1/dy.^2;
den= -2*(B+C); 

% Normalize Multipliers%
B = B/den; 
C = C/den; 
F = F/den; 
den = 1; 
error=10; 
error_iterations=0;
% check for diagonal dominance of elements 
abs(den) >= abs(2*B+2*C)
while error>10^-10; 
    W=U; 

for i = 2:N+1; 
     
    % Left boundary 
    U(i,1) = den*(  F(i,1) - (2*B)*U(i,2) - C*U(i-1,1) - C*U(i+1,1) );
    
    % Right Boundary 
    U(i,N+2) = den*(  F(i,N+2) - (2*B)*U(i,N+1) - C*U(i-1,N+2) - C*U(i+1,N+2) ); 
end 

%% Gauss-Siedel iterating the general U equation%

for i = 2:N+1; 
    for j = 2:M+1; 
        U(i,j) = den*(  F(i,j) - C*U(i+1,j) - C*U(i-1,j)- B*U(i,j+1) - B*U(i,j-1) ); 
    end 
end 
error=abs(max(max(((W-U)./W)))); 
error_iterations=error_iterations+1;
end 
toc 
error_iterations
figure 
subplot(1,2,1),surf(U),xlabel('x axis'),ylabel('y axis'),title('F=cos(x)sin(y)');

subplot(1,2,2),contour(U),xlabel('x axis'),ylabel('y axis'),title('F=cos(x)sin(y)');
