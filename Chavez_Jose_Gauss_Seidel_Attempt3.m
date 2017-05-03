%%%%%% Code Method 1 %%%%%%%% 
% Jose Chavez  1161146 
clear all; clc; 
tic 
%% Parameters 
ax = -pi; 
ay = -pi; 
bx = pi; 
by = pi; 
% Define the number of points on the interior (this does not include the 
% exterior boundary points) 
N=input('Value of X Intenal Nodes='); 
M=input('Value of Y Internal Nodes='); 
Me=M+2; 
Ne=N+2; 
% this generates the x and y values that will be used to calculate 
xvalues = linspace(-pi,pi,Ne); 
yvalues = linspace(-pi,pi,Me); 
%% 
%U matrix (guess) 
U = ones(Ne,Me); 
%solving for right hand side with F equation 
for i=1:length(xvalues); 
    for j=1:length(yvalues); 
F(i,j) = cos ( (0.5*pi)* (2*((xvalues(i)-ax) / (bx - ax))+1 )).*sin( pi*((yvalues(j)-ay) / (by -ay))); 
    end 
end 

%% Boundary Conditions for "top" and "bottom" 

% Bottom boundary values 
phi_ab = ((xvalues - ax).^2 ) .* sin( (pi *(xvalues - ax)) / (2*(bx-ax)) ) ; 

% Top boundary values 
psy_ab = cos (pi*(xvalues-ax)).*cosh(bx-xvalues); 

% place these known values in the solution grid 
U(:,1) = phi_ab; 
U(:,Me) = psy_ab; 
%% Left and Right Boundary points 
%   Using the given neumann condition yields special cases of the 
%   gauss-siedel iteration that can be used along entire "side" boundaries. 
%   F is already assumed to have been generated as well as U solution grid 
% Parameters that are used in the iterations. 
DX = 2*pi/(N+1); 
B = 1/DX.^2 
DY = 2*pi/(M+1); 
C = 1/DY.^2 
DEN = -2*(B+C) 

% normalize elements 
B = B/DEN; 
C = C/DEN; 
F = F/DEN; 
DEN = 1; 
error=10; 
error_iterations=0 
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
%% Main Sweep of Gauss-Siedel 

for i = 2:N+1; 
    for j = 2:M+1; 
        U(i,j) = DEN*(  F(i,j) - C*U(i+1,j) - C*U(i-1,j)- B*U(i,j+1) - B*U(i,j-1) ); 
    end 
end 
end 
error=abs(max(max(((W-U)./W)))); 
error_iterations=error_iterations+1 
end 
toc 
figure 
subplot(1,2,1),surf(U),xlabel('y axis'),ylabel('x axis');
subplot(1,2,2),contour(U),xlabel('y axis'),ylabel('x axis')
