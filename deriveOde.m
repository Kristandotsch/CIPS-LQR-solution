%pkg load symbolic
clc
clear all

% define symbolic variables
syms t;                 % time 
syms u;                 % control input
syms q1f(t);            % generalized koordinate which describes the cart position (symbolic function)
syms q2f(t);            % generalized koordinate which describes the pole angle (symbolic function)
syms g l m1 m2 d1;      % constants: gravitation, lenth of pole, mass of cart, mass of pendulum, viscous friction of the cart 

% mass matrices
M1=m1*eye(3); M1(3,3)=0;  % Moment of inertia = 0 -> cart doesn't turn
M2=m2*eye(3); M2(3,3)=0;  % Moment of inertia around cg = 0 -> concentrated mass at the end of the pole

% pose of the cart          
y1=[  q1f;... 
      0;...
      0     ];

% pose of the pole
y2=[  q1f-l*sin(q2f); ...
       l*cos(q2f);  ...
          q2f          ];

% derive poses
y1d=diff(y1,t);
y2d=diff(y2,t);

% calculate energies
T=0.5*(y1d.'*M1*y1d + y2d.'*M2*y2d);
gvec=[0; g; 0].';
U=gvec*M1*y1+gvec*M2*y2;

% calculate Lagrangian
L=T-U;
L=simplify(L);

% L_q: 
q=[q1f; q2f];
L_q=jacobian(L,q);   
L_q=simplify(L_q);

% L_qd:
qd=diff(q,t);
L_qd=jacobian(L,qd);
L_qd=simplify(L_qd);

% L_qd_t:
L_qd_t=diff(L_qd,t);
L_qd_t=simplify(L_qd_t);

% Euler-Lagrange formalism (note the external viscous friction force)
lhs=L_qd_t - L_q - [ u - d1 * diff(q1f,t) 0]; 
lhs=simplify(lhs)

% rename variables
syms q1dd q2dd q1d q2d q1 q2
lhs=subs(lhs,diff(diff(q1f,t),t),q1dd);
lhs=subs(lhs,diff(diff(q2f,t),t),q2dd);
lhs=subs(lhs,diff(q1f,t),q1d);
lhs=subs(lhs,diff(q2f,t),q2d);
lhs=subs(lhs,q1f(t),q1);
lhs=subs(lhs,q2f(t),q2);

% solve for the highest derrivatives
eq=lhs== 0;
Sol=solve(eq,q1dd,q2dd);
q1dd=simplify(Sol.q1dd);
p1ddCode=char(q1dd)
q2dd=simplify(Sol.q2dd);
p2ddCode=char(q2dd)

% state space representation: statevector = [x, x_dot, theta, theta_dot]
statevector=[q1; q1d; q2; q2d];
f=[q1d; q1dd; q2d; q2dd];         

% compute the jacobians
A=jacobian(f,statevector);                   
B=jacobian(f,u);

% evaluate jacobians at fixpoint = [0 0 0 0]
A_=subs(A,[q1,q1d, q2,q2d],[0,0,0,0])        
B_=subs(B,[q1,q1d, q2,q2d],[0,0,0,0])



                          
