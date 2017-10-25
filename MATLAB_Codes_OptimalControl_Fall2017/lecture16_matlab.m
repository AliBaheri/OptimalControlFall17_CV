%   MATLAB code for lecture 16's vehicle example
%   Written by Chris Vermillion

clear
clc

%   System parameters
m = 1000;
g = 9.8;
Crr = 0.01;
rho = 1.2;
Cd = 0.4;
A_ref = 5;

%   Initial conditions
x0 = 0;         %   Initial position (0)
v0 = 25;        %   Initial velocity (0)

%   Minumum and maximum allowable positions, velocities, and applied
%   forces, respectively
x_min = 0;      %   m
x_max = 1700;   %   m
v_min = 20;     %   m/s
v_max = 30;     %   m/s
u_min = 0;      %   N
u_max = 3000;   %   N

u_allowable = linspace(u_min,u_max,31);      %   Quantization of control forces (N)
x_allowable = linspace(x_min,x_max,41);     %   Quantization of positions (m)
v_allowable = linspace(v_min,v_max,21);     %   Quantization of velocities (m/s)

T_final = 60;   %   s
N = 30;         %   Number of stages
delta_T = T_final/N;    %   Time step

%   Populate a table of stage costs and destination states for quick reference
J_stage = zeros(length(x_allowable),length(v_allowable),length(u_allowable));
x_dest = zeros(length(x_allowable),length(v_allowable),length(u_allowable));
v_dest = zeros(length(x_allowable),length(v_allowable),length(u_allowable));

x_final = 1600;     %   Final state constraint
v_final = 25;       %   Minimum final velocity constraint

%   Terminal penalties
Kx_final_penalty = 100;
Kv_final_penalty = 1000000;

for i=1:length(x_allowable)
    for j=1:length(v_allowable)
        for k=1:length(u_allowable)
            v_dest(i,j,k) = v_allowable(j) + 1/m*(u_allowable(k) - Crr*m*g - 0.5*rho*Cd*A_ref*v_allowable(j)^2)*delta_T;
            x_dest(i,j,k) = x_allowable(i) + 0.5*(v_dest(i,j,k)+v_allowable(j))*delta_T;
            if x_dest(i,j,k)>=x_min && v_dest(i,j,k)>=v_min && v_dest(i,j,k)<=v_max        %   State constraint satisfied
                J_stage(i,j,k) = u_allowable(k)*v_allowable(j)*delta_T;
            else                            %   State constraint violated
                J_stage(i,j,k) = 10000000;
            end
            
            J_term(i,j,k) = u_allowable(k)*v_allowable(j)*delta_T + Kx_final_penalty*(max(0,x_final-x_dest(i,j,k)))^2 + Kv_final_penalty*(max(0,v_final-v_dest(i,j,k)))^2;
            %if x_dest(i,j,k)>=x_final && v_dest(i,j,k)>=v_final     %   Terminal state constraints satisfied
            %    J_term(i,j,k) = u_allowable(k)*v_allowable(j)*delta_T;
            %else                            %   Terminal state constraints violated
            %    J_term(i,j,k) = 10000000;
            %end
        end
    end
end

%   Compute the vector of optimal costs to go from step N-1, along with the
%   corresponding control values
stage = N-1
J_opt_togo = zeros(length(x_allowable),length(v_allowable));  %   Optimal cost to go from each state
u_opt_togo_matrix = zeros(length(x_allowable),length(v_allowable),N);    %   Optimal control sequence from each originating state, at each stage
for i=1:length(x_allowable)
    for j=1:length(v_allowable)
        [J_opt_togo(i,j),index] = min(J_term(i,j,:));
        u_opt_togo_matrix(i,j,N) = u_allowable(index);
    end
end

%   Begin the backward recursion through the rest of the stages
stage = N-2;
while stage >= 0
    stage
    J_opt_togo_prev = J_opt_togo;
    u_opt_togo_matrix_prev = u_opt_togo_matrix;
    for i=1:length(x_allowable)
        for j=1:length(v_allowable)
            for k=1:length(u_allowable)
                %   Calculate the intermediate (successor) state resulting from
                %   application of each candidate control input
                v_inter = v_allowable(j) + 1/m*(u_allowable(k) - Crr*m*g - 0.5*rho*Cd*A_ref*v_allowable(j)^2)*delta_T;
                x_inter = x_allowable(i) + 0.5*(v_allowable(j)+v_inter)*delta_T;
                %   In general, the successor state calculated above will NOT
                %   be one of the quantized state values, so we will need to
                %   interpolate
                if x_inter>=x_min && v_inter<=v_max && v_inter >=v_min
                    J_opt_togo_current = interp2(v_allowable,x_allowable,J_opt_togo_prev,v_inter,x_inter);
                    %   Calculate the optimal cost to go through the above
                    %   successor state
                    J_togo(k) = J_stage(i,j,k) + J_opt_togo_current;
                else
                    J_togo(k) = 1000000000;
                end
            end
            [J_opt_togo(i,j),control_index] = min(J_togo);
            u_opt_togo_matrix(i,j,stage+1) = u_allowable(control_index);
            %   Need to append the optimal control sequence with the optimal
            %   control sequence to go from the the intermediate state
            %   calculated above. Because this intermediate state will, in
            %   general, not correspond to one of the quantized states, we must
            %   interpolate the control sequence as well
            v_inter = v_allowable(j) + 1/m*(u_allowable(control_index) - Crr*m*g - 0.5*rho*Cd*A_ref*v_allowable(j)^2)*delta_T;
            v_inter = median([v_min v_inter v_max]);
            x_inter = x_allowable(i) + 0.5*(v_allowable(j)+v_inter)*delta_T;
            x_inter = median([x_min x_inter x_max]);
            for n=stage+2:N
                u_opt_togo_matrix(i,j,n) = interp2(v_allowable,x_allowable,u_opt_togo_matrix_prev(:,:,n),v_inter,x_inter);
            end
        end
    end
    stage = stage - 1;
end
%   Note - We have computed the optimal control input sequence for all
%   possible initial conditions. The example asks specifically about x0 =
%   1, so we will plot the optimal control sequence for that one in
%   particular
relevant_x0 = find(x_allowable == x0);
relevant_v0 = find(v_allowable == v0);

for i=1:N
    u_opt_vector(i) = u_opt_togo_matrix(relevant_x0,relevant_v0,i);
end

time_vector = linspace(0,60,N+1);

figure(1)
stairs(time_vector,[u_opt_vector u_opt_vector(length(u_opt_vector))]);
xlabel('Time (s)','fontsize',12);
ylabel('Control signal (N)','fontsize',12);

%   Report the optimal cost for the relevant IC
J_best_DP = J_opt_togo(relevant_x0,relevant_v0)

x_opt_vector(1) = x0;
v_opt_vector(1) = v0;
%   Reconstruct the state evolution under the optimized control sequence
for i=1:length(u_opt_vector)
    v_opt_vector(i+1) = v_opt_vector(i) + 1/m*(u_opt_vector(i) - Crr*m*g - 0.5*rho*Cd*A_ref*v_opt_vector(i)^2)*delta_T;
    x_opt_vector(i+1) = x_opt_vector(i) + 0.5*(v_opt_vector(i)+v_opt_vector(i+1))*delta_T;
end

figure(2)
stairs(time_vector,x_opt_vector);
xlabel('Time (s)','fontsize',12);
ylabel('Position (m)','fontsize',12);

figure(3)
stairs(time_vector,v_opt_vector);
xlabel('Time (s)','fontsize',12);
ylabel('Velocity (m/s)','fontsize',12);

u_opt_matrix = u_opt_togo_matrix;
