%   Lecture 18 m-file
%   By Chris Vermillion

clear

%   Define the radius and center of the obstacle to be avoided
R = 0.25;
xc = 0.5;
yc = 0;

%   Target and initial velocity parameters
x_target = 1;
y_target = 0;
v0 = 1;

%   Initial position and orientation
x_init = 0;
y_init = 0;
psi_init = 0.1;

%   Simulation parameters
end_time = 10;
T_MPC = 0.1;

sim('lecture18_simulink');

figure(1)
subplot(2,1,1)
hold on
plot(time_vec,x);
grid
xlabel('Time (s)','fontsize',12);
ylabel('x position (m)','fontsize',12);
subplot(2,1,2)
hold on
plot(time_vec,y);
grid
xlabel('Time (s)','fontsize',12);
ylabel('y position (m)','fontsize',12);

figure(2)
plot(x,y);
hold on
%   Draw a circle
param = linspace(0,2*pi,200);
for i=1:length(param)
    circle_x(i) = xc + R*cos(param(i));
    circle_y(i) = yc + R*sin(param(i));
end
plot(circle_x,circle_y,'r--');
legend('Path','Obstacle');
xlabel('x (m)','fontsize',12);
ylabel('y (m)','fontsize',12);
