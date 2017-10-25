%   MATLAB code for lecture 12
%   Written by Chris Vermillion

clear

%   Ex. 7.10 - One iteration of SQP
u0 = [2 .5]';       %   Initial guess
g0 = 1 - u0(1) - u0(2);     %   Initial value of the inequality constraint
h0 = 1 - u0(1)*u0(2);       %   Initial value of the equality constraint

%   Calculate the gradient and hessian of the objective function
del_J = [u0(2)^-2-u0(2)*u0(1)^-2 -2*u0(1)*u0(2)^-3+u0(1)^-1];
hessian = [2*u0(2)*u0(1)^-3 -2*u0(2)^-3-u0(1)^-2;
    -2*u0(2)^-3-u0(1)^-2 6*u0(1)*u0(2)^-4];

%   Calculate the gradients of the constraint functions
del_g = -[1 1];
del_h = -[u0(2) u0(1)];

iter_one = quadprog(hessian,del_J',del_g,-g0,del_h,-h0)

%   Example 2 from notes
u0 = [1 1]';        %   Initial guess
u_matrix = zeros(2,4);      %   This will contain guesses at each iteration
u_matrix(:,1) = u0;

for i=1:3
    del_J = [-exp(-u_matrix(1,i)) 2*(u_matrix(2,i)-2)];
    hessian = [exp(-u_matrix(1,i)) 0; 0 2];
    del_g = [u_matrix(2,i) u_matrix(1,i)];
    g = u_matrix(1,i)*u_matrix(2,i);
    
    delta_u = quadprog(hessian,del_J',del_g,1-g);
    
    u_new = u_matrix(:,i) + delta_u;
    u_matrix(:,i+1) = u_new;
end
u_matrix
