%   Example problems for lecture 7

clear

%   Example 1
clear J
clear J_iter
clear u

%   Set up a non-adaptive gradient descent optimization
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
figure(1)
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
figure(2)
plot(J_iter);
xlabel('Iteration','fontsize',12);
ylabel('J','fontsize',12);
title('Optimal control example - Gradient descent','fontsize',12);

%   Set up an adaptive gradient descent optimization, using a line search, for the same problem
clear J
clear J_iter
clear u
u(:,1) = [-1 -1 -1]';
J_iter(1) = 100 + u(1,1)^2 + (10+u(1,1))^2 + u(2,1)^2 + (10+u(1,1)+u(2,1))^2 + u(3,1)^2;
delta_u_convergence = .001;
delta_u = 1;
alpha_min = 0;
alpha_max = .5;
alpha_vec = linspace(alpha_min,alpha_max,100);
k = 1;
counter = 0;

%   Perform the gradient descent iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient(1) = 2*u(1,k) + 2*(10+u(1,k)) + 2*(10+u(1,k)+u(2,k));
    gradient(2) = 2*u(2,k) + 2*(10+u(1,k)+u(2,k));
    gradient(3) = 2*u(3,k);
    J_best = J_iter(k);
    u_best = u(:,k);
    for i=1:length(alpha_vec)
        u_prelim = u(:,k) - alpha_vec(i)*gradient';
        J_prelim = 100 + u_prelim(1)^2 + (10+u_prelim(1))^2 + u_prelim(2)^2 + (10+u_prelim(1)+u_prelim(2))^2 + u_prelim(3)^2;
        if J_prelim < J_best
            u_best = u_prelim;
            J_best = J_prelim;
        end
    end
    u(:,k+1) = u_best;
    J_iter(k+1) = 100 + u(1,k+1)^2 + (10+u(1,k+1))^2 + u(2,k+1)^2 + (10+u(1,k+1)+u(2,k+1))^2 + u(3,k+1)^2;
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the convergence of the elements of the control signal
figure(3)
plot(u(1,:));
hold on
plot(u(2,:),'r');
hold on
plot(u(3,:),'m');
xlabel('Iteration','fontsize',12);
ylabel('Control signals','fontsize',12);
legend('u(1)','u(2)','u(3)');
title('Optimal control example - Modified gradient descent with line search','fontsize',12);

%   Plot the evolution of the cost function
figure(4)
plot(J_iter);
xlabel('Iteration','fontsize',12);
ylabel('J','fontsize',12);
title('Optimal control example - Modified gradient descent with line search','fontsize',12);

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

figure(5)
mesh(u1,u2,J)
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Actual objective function vs. design variables','fontsize',12);

%   Set up a (pure) Netwon's method optimization
clear J
clear J_iter
clear u
u(:,1) = [2 -2]';
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
title('Pure Newton-based optimization','fontsize',12);

%   Set up a modified Netwon's method optimization with line search
clear J
clear J_iter
clear u
u(:,1) = [2 -2]';
J_iter(1) = u(1,1) + u(2,1) + u(1,1)*exp(-u(2,1)) + u(2,1)^2*exp(-u(1,1));
delta_u_convergence = .001;
delta_u = 1;
alpha_min = 0;
alpha_max = 1;
k = 1;
counter = 0;

alpha_vec = linspace(alpha_min,alpha_max,100);

%   Perform the Newton iterations
while delta_u >= delta_u_convergence && counter <= 100
    gradient = [1-exp(-u(2,k))-u(2,k)^2*exp(-u(1,k)) 1-u(1,k)*exp(-u(2,k))+2*u(2,k)*exp(-u(1,k))];
    hessian = [-u(2,k)^2*exp(-u(1,k)) -exp(-u(2,k))-2*u(2,k)*exp(-u(1,k)); ...
        -exp(-u(2,k))-2*u(2,k)*exp(-u(1,k)) u(1,k)*exp(-u(2,k))+2*exp(-u(1,k))];
    newton_direction = -inv(hessian)*gradient';
    u_best = u(:,k);
    J_best = J_iter(k);
    for i= 1:length(alpha_vec)
        u_prelim = u(:,k) + alpha_vec(i)*newton_direction;
        J_prelim = u_prelim(1) + u_prelim(2) + u_prelim(1)*exp(-u_prelim(2)) + u_prelim(2)^2*exp(-u_prelim(1));
        if J_prelim < J_best
            u_best = u_prelim;
            J_best = J_prelim;
        end
    end
    u(:,k+1) = u_best;
    J_iter(k+1) = u(1,k+1) + u(2,k+1) + u(1,k+1)*exp(-u(2,k+1)) + u(2,k+1)^2*exp(-u(1,k+1));
    delta_u = norm(u(:,k+1)-u(:,k));
    counter = counter + 1;
    k = k + 1;
end

%   Plot the path of the Newton's method-based optimization
figure(7)
plot3(u(1,:)',u(2,:)',J_iter);
xlabel('u1','fontsize',12);
ylabel('u2','fontsize',12);
zlabel('J','fontsize',12);
title('Modified Newton descent with line search','fontsize',12);
