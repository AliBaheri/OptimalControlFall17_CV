%   MATLAB code for lecture 11
%   Written by Chris Vermillion

clear

%   State space representation - note that C doesn't matter
A = [.99 .1; -.2 .9];
B = [0 .1]';

%   Initial condition and horizon length
x0 = [5 0]';
N = 50;

%   LQR weighting matrices
Q_bar = [1 0; 0 0];
r_bar = 1;

%   Compute the entries to the lifted system matrix
H = zeros(2*N,N);         %   Initialize H
for i=1:N
    for j=1:i
        H(2*(i-1)+1:2*i,j) = A^(i-j)*B;
    end
end

%   Compute Q
Q = zeros(N);         %   Initialize Q
for i=1:N
    Q = Q + H(2*(i-1)+1:2*i,:)'*Q_bar*H(2*(i-1)+1:2*i,:);
end
Q = Q + r_bar*eye(N);

%   Compute r
r = zeros(1,N);      %   Initialize r
for i=1:N
    r = r + 2*x0'*(A^i)'*Q_bar*H(2*(i-1)+1:2*i,:);
end

%   Identify constraints
A1 = [eye(N);
    -eye(N)];
b1 = 0.5*ones(2*N,1);

%   Use quadprog to determine the optimal control sequence
[u_opt,J_opt] = quadprog(2*Q,r,A1,b1,[],[]);

figure(1)
stairs(u_opt);
xlabel('Step','fontsize',12);
ylabel('Control signal','fontsize',12);

%   Design a linear feedback controller for comparison
K = dlqr(A,B,Q_bar,r_bar);

x = x0;
for i=1:N
    u(i) = -K*x;
    x = A*x + B*u(i);
end

figure(2)
stairs(u);
xlabel('Step','fontsize',12);
ylabel('Linear feedback control signal','fontsize',12);
