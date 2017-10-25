%   MATLAB code for lecture 13
%   Written by Chris Vermillion

clear

%   Example 1 - The example that yielded problems when the
%   Hessian of the objective function was used in the QP subproblem

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');   %   Last entry sets the algorithm to SQP
u0 = [1.2 .2]';         %   Initial guess
obj_fun = @ex1_lec13_objective;     %   Name of objective function
nl_constraint = @ex1_lec13_constraint;      %   Name of constraint function

u_opt_ex1 = fmincon(obj_fun,u0,[],[],[],[],[],[],nl_constraint,options)

%   Example 2 - Second example of lecture 12, revisited

u0 = [1 1]';            %   Initial guess
obj_fun = @ex2_lec13_objective;     %   Name of objective function
nl_constraint = @ex2_lec13_constraint;      %   Name of constraint function

u_opt_ex2 = fmincon(obj_fun,u0,[],[],[],[],[],[],nl_constraint,options)
