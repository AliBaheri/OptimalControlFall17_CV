%   MATLAB code for lecture 10
%   Written by Chris Vermillion

clear

%   Example 1

%   Encode constraints
A_constraint = [1 2; 1 -1; -2 0; 0 -2];
b_constraint = [8 1.5 -1 -1]';

%   Encode the gradient of the objective function
k = [-2 1]';

%   Initially consider the solution at the intersection of constraints 3
%   and 4
A_small = A_constraint(3:4,:);
b_small = [b_constraint(3) b_constraint(4)]';
u_init = inv(A_small)*b_small;
J_init = -2*u_init(1) + u_init(2);

convergence_flag = 0;
active_constraint_set = [3 4];
j = 1;
for i=1:length(b_constraint)
    if active_constraint_set(1)~=i && active_constraint_set(2)~=i
        inactive_constraint_set(j) = i;
        j = j + 1;
    end
end

while convergence_flag == 0
    %   Determine which of the active constraints to search along (this is
    %   equivalent to identifying the pivot column when doing LP by hand
    %   (tabular form). Sometimes this is done based on marginal utility;
    %   however, it is sufficient to choose *any* constraint direction
    %   along which the cost function is decreasing
    solution_found = 0;
    i = 1;
    while solution_found == 0 && i <= length(active_constraint_set)
        % By choosing an active constraint direction out of two options, we
        % are deleting one of the constraints. We must choose another
        % constraint to add, which is geometrically equivalent to finding
        % which vertix we want to evaluate the objective function at next.
        % We always choose the vertex closest to where we are, so long as
        % this leads to a decreasing cost function value. This step is
        % equivalent to choosing the pivot row when doing LP by hand
        % (tabular form).
        active_constraint = i;
        for i=1:length(inactive_constraint_set)
            A_small = [A_constraint(active_constraint,:); A_constraint(inactive_constraint_set(i),:)];
            b_small = [b_constraint(active_constraint,:); b_constraint(inactive_constraint_set(i),:)];
            u_candidate(:,i) = inv(A_small)*b_small;
            if k'*(u_candidate(:,i)-u_init) < 0       % Direction of decreasing cost
                dist(i) = norm(u_candidate(:,i) - u_init);
                solution_found = 1;
            else
                dist(i) = 1000+i;
            end
        end
        if solution_found == 1              %   Change the new active constraint set and get ready for another iteration
            [min_dist new_constraint] = min(dist);
            u_init = u_candidate(:,new_constraint);
            J_init = -2*u_init(1) + u_init(2);
            
            active_constraint_set = [active_constraint new_constraint];
            j = 1;
            for i=1:length(b_constraint)
                if active_constraint_set(1)~=i && active_constraint_set(2)~=i
                    inactive_constraint_set(j) = i;
                    j = j + 1;
                end
            end
        end
    end
    if solution_found ~= 0  %   This means that the cost function did not increase along any feasible constraint boundary, so u_init must be the minimizer
        convergence_flag = 1;
        u_opt = u_init;
        J_opt = J_init;
    end
end

u_opt
J_opt

%   Now let's do the LP the easy way...using the built-in MATLAB function!

%   Choose the dual simplex algorithm (dual because we're doing a
%   minimization)
options = optimoptions('linprog','Algorithm','dual-simplex');

[u_opt,J_opt] = linprog(k,A_constraint,b_constraint,[],[],[],[],options);

u_opt
J_opt
