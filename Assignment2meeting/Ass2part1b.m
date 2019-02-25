%% Part 1.2: Electrostatic Potential in 2D region
% The second stage involved modifying the boundary contitions of the
% region. The boundary conditions were set such that $V = V_0$ at $x = 0, x
% = L$ and at $y = 0, and y = W$, $V = 0$. The plot is generated and
% compared to the analytical solution. 
close all
clear

L = 3;
W = 2;
V_0 = 1;

dx = 0.05;
dy = 0.05;

nx = L/dx;
ny = W/dy;

c1 = 1/(dx^2); 
c2 = 1/(dy^2);
c3 = -2*(1/dx^2 + 1/dy^2);

G = zeros(nx*ny,nx*ny);

for i = 2:nx-1
    for j = 2:ny-1
        
        n = i + (j - 1)*nx;
        nym = i + (j - 2)*nx;
        nyp = i + j*nx;
        nxm = (i - 1) + (j - 1)*nx;
        nxp = (i + 1) + (j - 1)*nx;
        
        G(n,n) = c3;
        G(n,nxm) = c1;
        G(n,nxp) = c1;
        G(n,nym) = c2;
        G(n,nyp) = c2;
    end
end

F = zeros(nx*ny,1);


%% 
% Given the changes to the boundary conditions, the matrices are filled
% using the iterative approach. 

for j = 1:ny
    n = 1 + (j - 1)*nx;
    
    G(n, n) = 1;
    F(n) = V_0;
    
    n1 = nx + (j-1)*nx;
    G(n1, n1) = 1;
    F(n1) = 1;  
end

for i = 1:nx
    
    n = i;
    G(n,n) = 1;
    
    n1 = i + (ny - 1)*nx;
    G(n1, n1) = 1;
    
    
end

F(1) = 0;
F(1 + (ny - 1)*nx) = 0;
F(nx) = 0;
F(nx + (ny - 1)*nx) = 0;

%Finding solution and reshaping the transpose
V = G\F;
solution = reshape(V,[nx, ny])';

x = linspace(0, L, nx);
y = linspace(0, W, ny);

figure(1)
surf(x, y, solution)
xlabel('x')
ylabel('y')
title('Finite Differences Solution')

%% 
% The analytical solution was then determined and compared to the solution
% generated for the Finite Difference Solution. The analyitcal solution was
% calculated by using an infinite series provided in the assignment. 

Asol = zeros(ny, nx);
x_anyl = repmat(linspace(-L/2,L/2,nx),ny,1);
y_anyl = repmat(linspace(0,W,ny),nx,1)';
count = 100;

for i=1:count
    count = 2*i - 1;
    Asol = Asol + 1./count.*cosh(count.*pi.*x_anyl./W)./cosh(count.*pi.*(L./2)./W).*sin(count.*pi.*y_anyl./W);
end

Asol = Asol.*4.*V_0./pi;

%%
% The analytical solution can be seen in the figure below. It was noted
% that these solutions are near identical. The different between the two
% stems from the time taken to reach a reasonable solution. 


figure(2);
surf(linspace(0,L,nx),linspace(0,W,ny),Asol);
xlabel('x');
ylabel('y');
title('Analytical Solution');

%% 
% Generating a movie to compare the two solutions, analytical and using
% Finite Difference, displays the difference between the two. The variable
% dependance of the analytical solution results in the plot taking longer
% to take shape. The analytical solution can be stopped when it is within a
% small error margin of the Finite Difference solution. For more complex
% solutions, numerical method is advantageous as it is easier to apply
% compared to the analytical approach. In simple systems such as the
% rectangular region in this assignment, the analytical solution is valid. 


        