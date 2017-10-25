%   MATLAB code for lecture 14
%   Written by Chris Vermillion

clear
clc

%   Given the system description in the example problem, compute the
%   lifted system representation
A = 2;
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
N_allowable = length(u_allowable);
J_best = 10000; %   Initialize the best objective function so far to something ridiculously big

%   Search through all allowable control sequences
for i=1:N_allowable
    for j=1:N_allowable
        for k=1:N_allowable
            for m=1:N_allowable
                for n=1:N_allowable
                    u_sequence = [u_allowable(i) u_allowable(j) u_allowable(k) u_allowable(m) u_allowable(n)]';
                    x_sequence = H*u_sequence + M*x0;
                    J = 5*sum(x_sequence.^2) + sum(u_sequence.^2);
                    if J < J_best
                        J_best = J;
                        i_best = i;
                        j_best = j;
                        k_best = k;
                        m_best = m;
                        n_best = n;
                    end
                end
            end
        end
    end
end
u_best = [u_allowable(i_best) u_allowable(j_best) u_allowable(k_best) u_allowable(m_best) u_allowable(n_best)]';
x_best = H*u_best + M*x0;

figure(1)
stairs(u_best);

figure(2)
stairs(x_best);

%   Report the optimal cost
J_best

%%
%   Repeat the example using dynamic programming with backward recursion
x_allowable = [-1 -0.5 0 0.5 1];

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
state_index = ones(length(u_allowable),1);
while stage >= 0
    u_opt_togo_matrix_prev = u_opt_togo_matrix;
    for i=1:length(x_allowable)
        for j=1:length(u_allowable)
            %   Calculate the intermediate (successor) state resulting from
            %   application of each candidate control input
            x_inter = A*x_allowable(i) + B*u_allowable(j);
            if abs(x_inter) <= 1
                state_index(j) = find(x_allowable == x_inter);
                %   Calculate the optimal cost to go through the above
                %   successor state
                J_togo(j) = J_stage(i,j) + J_opt_togo(state_index(j));
            else
                J_togo(j) = 10000;
            end
        end
        [J_opt_togo(i),control_index] = min(J_togo);
        u_opt_togo_matrix(i,stage+1) = u_allowable(control_index);
        %   Need to append the optimal control sequence with the optimal
        %   control sequence to go from the the intermediate state
        %   calculated above
        u_opt_togo_matrix(i,stage+2:N) = u_opt_togo_matrix_prev(state_index(control_index),stage+2:N);
    end
    stage = stage - 1;
end
%   Note - We have computed the optimal control input sequence for all
%   possible initial conditions. The example asks specifically about x0 =
%   1, so we will plot the optimal control sequence for that one in
%   particular
relevant_IC = find(x_allowable == x0);
u_opt_vector = u_opt_togo_matrix(relevant_IC,:);

figure(3)
stairs(u_opt_vector);

%   Report the optimal cost for the relevant IC
J_best_DP = J_opt_togo(relevant_IC)

%   Reconstruct the state evolution under the optimized control sequence
x_opt_vector = H*u_opt_vector' + M*x0;

figure(4)
stairs(x_opt_vector);

%   Although the plots in figures 3 and 4 are given for x0 = 1, we have in
%   fact computed the optimal control sequences for all candidate initial
%   conditions (-1, -0.5, 0, 0.5, and 1). These are given in the matrix
%   u_opt_togo_matrix
u_opt_matrix = u_opt_togo_matrix
