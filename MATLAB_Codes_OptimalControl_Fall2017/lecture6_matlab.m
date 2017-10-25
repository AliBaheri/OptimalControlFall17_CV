%   Example problems for lecture 6

clear

%   Example 1

%   First, grid a 2D domain and plot the cost function to visualize the
%   response
u1 = linspace(-2,2,40);
u2 = linspace(-2,2,40);

for i=1:length(u1)
    for j=1:length(u2)
        J(i,j) = 4*u1(i)^2 + 3*u1(i)*u2(j) + u2(j)^2;
    end
end

figure(1)
hold on
mesh(u1,u2,J)
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 1 - Actual objective function vs. design variables','fontsize',12);

%   Set up a gradient descent optimization
u(:,1) = [1 1]';
J_iter(1) = 4*u(1,1)^2 + 3*u(1,1)*u(2,1) + u(2,1)^2;
delta_u_convergence = .001;
delta_u = 1;
alpha = .1;
k = 1;
counter = 0;

%   Perform the gradient descent iteration
while delta_u >= delta_u_convergence && counter <= 100
    gradient = [8*u(1,k)+3*u(2,k) 3*u(1,k)+2*u(2,k)];
    u(:,k+1) = u(:,k) - alpha*gradient';
    J_iter(k+1) = 4*u(1,k+1)^2 + 3*u(1,k+1)*u(2,k+1) + u(2,k+1)^2;
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the path of the gradient descent optimization
figure(2)
hold on
plot3(u(1,:)',u(2,:)',J_iter);
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 1 - Gradient descent optimization','fontsize',12);

%   Set up a (pure) Netwon's method optimization for the same problem
clear J
clear J_iter
clear u

u(:,1) = [1 1]';
J_iter(1) = 4*u(1,1)^2 + 3*u(1,1)*u(2,1) + u(2,1)^2;
delta_u_convergence = .001;
delta_u = 1;
k = 1;
counter = 0;

%   Perform the Newton iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient = [8*u(1,k)+3*u(2,k) 3*u(1,k)+2*u(2,k)];
    hessian = [8 3; 3 2];
    u(:,k+1) = u(:,k) - inv(hessian)*gradient';
    J_iter(k+1) = 4*u(1,k+1)^2 + 3*u(1,k+1)*u(2,k+1) + u(2,k+1)^2;
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the path of the Newton-based optimization
figure(3)
plot3(u(1,:)',u(2,:)',J_iter);
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 1 - Pure Newton-based optimization','fontsize',12);

%   Example 2
clear J
clear J_iter
clear u

u1 = linspace(-3,3,40);
u2 = linspace(-3,3,40);

for i=1:length(u1)
    for j=1:length(u2)
        J(i,j) = u1(i) + u2(j) + u1(i)*exp(-u2(j)) + u2(j)^2*exp(-u1(i));
    end
end

figure(4)
mesh(u1,u2,J)
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 2 - Actual objective function vs. design variables','fontsize',12);

%   Set up a gradient descent optimization
u(:,1) = [-2 2]';
J_iter(1) = u(1,1) + u(2,1) + u(1,1)*exp(-u(2,1)) + u(2,1)^2*exp(-u(1,1));
delta_u_convergence = .001;
delta_u = 1;
alpha = .01;
k = 1;
counter = 0;

%   Perform the gradient descent iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient = [1-exp(-u(2,k))-u(2,k)^2*exp(-u(1,k)) 1-u(1,k)*exp(-u(2,k))+2*u(2,k)*exp(-u(1,k))];
    u(:,k+1) = u(:,k) - alpha*gradient';
    J_iter(k+1) = u(1,k+1) + u(2,k+1) + u(1,k+1)*exp(-u(2,k+1)) + u(2,k+1)^2*exp(-u(1,k+1));
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the path of the gradient descent optimization
figure(5)
plot3(u(1,:)',u(2,:)',J_iter);
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 2 - Gradient descent optimization','fontsize',12);

%   Set up a (pure) Netwon's method optimization for the same problem
clear J
clear J_iter
clear u
u(:,1) = [-2 2]';
J_iter(1) = u(1,1) + u(2,1) + u(1,1)*exp(-u(2,1)) + u(2,1)^2*exp(-u(1,1));
delta_u_convergence = .001;
delta_u = 1;
k = 1;
counter = 0;

%   Perform the Newton iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient = [1-exp(-u(2,k))-u(2,k)^2*exp(-u(1,k)) 1-u(1,k)*exp(-u(2,k))+2*u(2,k)*exp(-u(1,k))];
    hessian = [-u(2,k)^2*exp(-u(1,k)) -exp(-u(2,k))-2*u(2,k)*exp(-u(1,k)); ...
        -exp(-u(2,k))-2*u(2,k)*exp(-u(1,k)) u(1,k)*exp(-u(2,k))+2*exp(-u(1,k))];
    u(:,k+1) = u(:,k) - inv(hessian)*gradient';
    J_iter(k+1) = u(1,k+1) + u(2,k+1) + u(1,k+1)*exp(-u(2,k+1)) + u(2,k+1)^2*exp(-u(1,k+1));
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the path of the Newton's method-based optimization
figure(6)
plot3(u(1,:)',u(2,:)',J_iter);
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Example 2 - Pure Newton-based optimization','fontsize',12);

%   Example 3
clear J
clear J_iter
clear u

%   Set up a gradient descent optimization
u(:,1) = [-1 -1 -1]';
J_iter(1) = 100 + u(1,1)^2 + (10+u(1,1))^2 + u(2,1)^2 + (10+u(1,1)+u(2,1))^2 + u(3,1)^2;
delta_u_convergence = .001;
delta_u = 1;
alpha = .1;
k = 1;
counter = 0;

%   Perform the gradient descent iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient(1) = 2*u(1,k) + 2*(10+u(1,k)) + 2*(10+u(1,k)+u(2,k));
    gradient(2) = 2*u(2,k) + 2*(10+u(1,k)+u(2,k));
    gradient(3) = 2*u(3,k);
    u(:,k+1) = u(:,k) - alpha*gradient';
    J_iter(k+1) = 100 + u(1,k+1)^2 + (10+u(1,k+1))^2 + u(2,k+1)^2 + (10+u(1,k+1)+u(2,k+1))^2 + u(3,k+1)^2;
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the convergence of the elements of the control signal
figure(7)
plot(u(1,:));
hold on
plot(u(2,:),'r');
hold on
plot(u(3,:),'m');
xlabel('Iteration','fontsize',12);
ylabel('Control signals','fontsize',12);
legend('u(1)','u(2)','u(3)');
title('Optimal control example - Gradient descent','fontsize',12);

%   Plot the evolution of the cost function
figure(8)
plot(J_iter);
xlabel('Iteration','fontsize',12);
ylabel('J','fontsize',12);
title('Optimal control example - Gradient descent','fontsize',12);

%   Set up a (pure) Netwon's method optimization for the same problem
clear J
clear J_iter
clear u
u(:,1) = [-1 -1 -1]';
J_iter(1) = 100 + u(1,1)^2 + (10+u(1,1))^2 + u(2,1)^2 + (10+u(1,1)+u(2,1))^2 + u(3,1)^2;
delta_u_convergence = .001;
delta_u = 1;
k = 1;
counter = 0;

%   Perform the Newton iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient(1) = 2*u(1,k) + 2*(10+u(1,k)) + 2*(10+u(1,k)+u(2,k));
    gradient(2) = 2*u(2,k) + 2*(10+u(1,k)+u(2,k));
    gradient(3) = 2*u(3,k);
    hessian = [6 2 0; 2 4 0; 0 0 2];
    u(:,k+1) = u(:,k) - inv(hessian)*gradient';
    J_iter(k+1) = 100 + u(1,k+1)^2 + (10+u(1,k+1))^2 + u(2,k+1)^2 + (10+u(1,k+1)+u(2,k+1))^2 + u(3,k+1)^2;
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the convergence of the elements of the control signal
figure(9)
plot(u(1,:));
hold on
plot(u(2,:),'r');
hold on
plot(u(3,:),'m');
xlabel('Iteration','fontsize',12);
ylabel('Control signals','fontsize',12);
legend('u(1)','u(2)','u(3)');
title('Optimal control example - Pure Newton-based optimization','fontsize',12);

%   Plot the evolution of the cost function
figure(10)
plot(J_iter);
xlabel('Iteration','fontsize',12);
ylabel('J','fontsize',12);
title('Optimal control example - Pure Newton-based optimization','fontsize',12);
