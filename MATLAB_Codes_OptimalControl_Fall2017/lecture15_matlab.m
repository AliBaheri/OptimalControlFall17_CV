%   MATLAB code for lecture 15
%   Written by Chris Vermillion

clear
clc

%   Given the system description in the example problem, compute the
%   lifted system representation
A = 1.5;
B = 0.5;
N = 5;
H = zeros(N,N);         %   Initialize H
M = zeros(N,1);
for i=1:N
    for j=1:i
        H(i,j) = A^(i-j)*B;
    end
    M(i) = A^i;
end
x0 = 1;         %   Initial condition

u_allowable = [-3 -2 -1 0 1 2 3];   %   Quantization of control values
x_allowable = [-1 -0.5 0 0.5 1];    %   Quantization of state values

%   Populate a table of stage costs and destination states for quick reference
J_stage = zeros(length(x_allowable),length(u_allowable));
x_dest = zeros(length(x_allowable),length(u_allowable));
for i=1:length(x_allowable)
    for j=1:length(u_allowable)
        x_dest(i,j) = A*x_allowable(i) + B*u_allowable(j);
        if abs(x_dest(i,j)) <= 1        %   State constraint satisfied
            J_stage(i,j) = 5*x_dest(i,j)^2 + u_allowable(j)^2;
        else                            %   State constraint violated
            J_stage(i,j) = 10000;
        end
    end
end

%   Compute the vector of optimal costs to go from step N-1, along with the
%   corresponding control values
J_opt_togo = zeros(length(x_allowable),1);  %   Optimal cost to go from each state
u_opt_togo_matrix = zeros(length(x_allowable),N);    %   Optimal control sequence from each originating state, at each stage
for i=1:length(x_allowable)
    [J_opt_togo(i),index] = min(J_stage(i,:));
    u_opt_togo_matrix(i,N) = u_allowable(index);
end

%   Begin the backward recursion through the rest of the stages
stage = N-2;
while stage >= 0
    J_opt_togo_prev = J_opt_togo;
    u_opt_togo_matrix_prev = u_opt_togo_matrix;
    for i=1:length(x_allowable)
        for j=1:length(u_allowable)
            %   Calculate the intermediate (successor) state resulting from
            %   application of each candidate control input
            x_inter = A*x_allowable(i) + B*u_allowable(j);
            %   In general, the successor state calculated above will NOT
            %   be one of the quantized state values, so we will need to
            %   interpolate
            if abs(x_inter) <= 1
                J_opt_togo_current = interp1(x_allowable,J_opt_togo_prev,x_inter);
                %   Calculate the optimal cost to go through the above
                %   successor state
                J_togo(j) = J_stage(i,j) + J_opt_togo_current;
            else
                J_togo(j) = 10000;
            end
        end
        [J_opt_togo(i),control_index] = min(J_togo);
        u_opt_togo_matrix(i,stage+1) = u_allowable(control_index);
        %   Need to append the optimal control sequence with the optimal
        %   control sequence to go from the the intermediate state
        %   calculated above. Because this intermediate state will, in
        %   general, not correspond to one of the quantized states, we must
        %   interpolate the control sequence as well
        x_inter = A*x_allowable(i) + B*u_allowable(control_index);
        for k=stage+2:N
            u_opt_togo_matrix(i,k) = interp1(x_allowable,u_opt_togo_matrix_prev(:,k),x_inter);
        end
    end
    stage = stage - 1;
end
%   Note - We have computed the optimal control input sequence for all
%   possible initial conditions. The example asks specifically about x0 =
%   1, so we will plot the optimal control sequence for that one in
%   particular
relevant_IC = find(x_allowable == x0);
u_opt_vector = u_opt_togo_matrix(relevant_IC,:);

figure(1)
stairs(u_opt_vector);
xlabel('Time step','fontsize',12);
ylabel('u','fontsize',12);

%   Report the optimal cost for the relevant IC
J_best_DP = J_opt_togo(relevant_IC)

%   Reconstruct the state evolution under the optimized control sequence
x_opt_vector = H*u_opt_vector' + M*x0;

figure(2)
stairs(x_opt_vector);
xlabel('Time step','fontsize',12);
ylabel('x','fontsize',12);

%   Although the plots in figures 3 and 4 are given for x0 = 1, we have in
%   fact computed the optimal control sequences for all candidate initial
%   conditions (-1, -0.5, 0, 0.5, and 1). These are given in the matrix
%   u_opt_togo_matrix
u_opt_matrix = u_opt_togo_matrix

%   Repeat the exercise using forward recursion

%   Compute the vector of optimal costs to arrive at step 1, along with the
%   corresponding control values
J_opt_toarrive = zeros(length(x_allowable),1);  %   Optimal cost to arrive at each state
u_opt_toarrive_matrix = zeros(length(x_allowable),N);    %   Optimal control sequence to each terminal state, at each stage
for i=1:length(x_allowable)
    %   Determine the required control input such that the originating
    %   state matches the initial condition (x0)
    u_req = 1/B*(x_allowable(i) - A*x0);
    orig_index = find(x_allowable == x0);
    if abs(u_req) <= 3
        J_opt_toarrive(i) = interp1(u_allowable,J_stage(orig_index,:),u_req);
    else
        J_opt_toarrive(i) = 10000;
    end
    u_opt_toarrive_matrix(i,1) = u_req;
end

%   Begin the forward recursion through the rest of the stages
stage = 2;
while stage <= N
    J_opt_toarrive_prev = J_opt_toarrive;
    u_opt_toarrive_matrix_prev = u_opt_toarrive_matrix;
    for i=1:length(x_allowable)
        for j=1:length(u_allowable)
            %   Calculate the previous state before the application of each
            %   candidate control input
            x_prev = inv(A)*(x_allowable(i) - B*u_allowable(j));
            %   In general, the previous state calculated above will NOT
            %   be one of the quantized state values, so we will need to
            %   interpolate
            if abs(x_prev) <= 1
                J_opt_toarrive_current = interp1(x_allowable,J_opt_toarrive_prev,x_prev);
                %   Calculate the optimal cost to go through the above
                %   successor state
                J_stage_current = 5*x_allowable(i)^2 + u_allowable(j)^2;
                J_toarrive(j) = J_stage_current + J_opt_toarrive_current;
            else
                J_toarrive(j) = 10000;
            end
        end
        [J_opt_toarrive(i),control_index] = min(J_toarrive);
        u_opt_toarrive_matrix(i,stage) = u_allowable(control_index);
        %   Need to append the optimal control sequence with the optimal
        %   control sequence to go from the the intermediate state
        %   calculated above. Because this intermediate state will, in
        %   general, not correspond to one of the quantized states, we must
        %   interpolate the control sequence as well
        x_prev = inv(A)*(x_allowable(i) - B*u_allowable(control_index));
        for k=1:stage-1
            u_opt_toarrive_matrix(i,k) = interp1(x_allowable,u_opt_toarrive_matrix_prev(:,k),x_prev);
        end
    end
    stage = stage + 1;
end
%   Note - We have computed the optimal control input sequence for all
%   possible terminal conditions. Let's look at the optimal control
%   sequence for x_final = 0 (recall that x0 = 1)
xf = 0;
relevant_TC = find(x_allowable == xf);
u_opt_vector = u_opt_toarrive_matrix(relevant_TC,:);

figure(3)
stairs(u_opt_vector);
xlabel('Time step','fontsize',12);
ylabel('u','fontsize',12);

%   Report the optimal cost for the relevant IC
J_best_DP = J_opt_toarrive(relevant_TC)

%   Reconstruct the state evolution under the optimized control sequence
x_opt_vector = H*u_opt_vector' + M*x0;

figure(4)
stairs(x_opt_vector);
xlabel('Time step','fontsize',12);
ylabel('x','fontsize',12);

%   Although the plots in figures 3 and 4 are given for xf = 0, we have in
%   fact computed the optimal control sequences for all candidate terminal
%   conditions (-1, -0.5, 0, 0.5, and 1). These are given in the matrix
%   u_opt_toarrive_matrix
u_opt_matrix = u_opt_toarrive_matrix
