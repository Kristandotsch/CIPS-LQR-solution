clear all 
close all 
clc
%pkg load control

% parameterization
m1 = 1;          % mass of the cart
m2 = 0.5;        % mass of the pole
l = 0.7;         % length of the pole
g = 9.81;        % earth gravity
d1 = 0.8;          % viscous friction on the cart 

% Jacobians:
A = [0 1 0 0;
    0 -d1/m1  m2*g/m1 0;
    0 0 0 1;
    0 -d1/(m1*l) (m1+m2)*g/(m1*l) 0];
    
B = [0; 1/m1; 0; 1/(m1*l)];

%compute eigenvalues
eig(A);

% check rank to see if the system is controllable
rank(ctrb(A,B));

% weight matrices for the optimization problem
Q = [100 0 0 0;   % x_ = [x x_d theta theta_d]
      0 1 0 0;
      0 0 1000 0;
      0 0 0 1];
R = 0.0001;

% solve algebraic Riccati equation to obtain gain vector 
%K = lqr(A,B,Q,R);
K = [-1000.00   -946.79   4353.27    803.53];

% create time vector
tspan = 0:0.01:15;

% initial condition: [x, x_d, theta, theta_d]
x0 = [0; 0; -pi; 2]; 
%x0 = [-1; 0; 0; 0.5];

% ode options
odeset('relTol', 1e-7, 'absTol', realmin);

% numericaly solve the ode
[t_out,x_out] = ode45( @(t,x)evaluateOde(t,x,m1,m2,l,g,d1,K), tspan, x0, odeset);

% draw the pendulum animation
drawPendulum(t_out, x_out, m1, m2, g, l);

