%% Correct implementation of the Newton-Euler Inverse Dynamics Algorithm
clear all;
close all;
clc;

theta = [0;pi/2];
theta_dot = [0.2;0.2];
% theta_double_dot = zeros(size(theta,1),1);
omega = [0;0;1];
Ftip = zeros(6,1);
g = [0;-9.8;0];

n = size(theta,1);
V0 = zeros(6,1);
Vdot_0 = zeros(6,1);
Vdot_0(1:3) = -g;
% This value is just to recheck and confirm the working of one of the
% functions
tau = zeros(size(theta));

m1 = 1; % Kg
a1 = 0.4; % m
lc1 = a1/2; % m

m2 = 1; % Kg
a2 = 0.4; % m
lc2 = a2/2; % m

% Rotation matrix and position vectors when all the joint variables are zero. 
R = eye(3);
p_00 = [0;0;0];
p_01 = [lc1;0;0];
p_02 = [a1 + lc1;0;0];
p_03 = [a1 + a2;0;0];

% Transformation matrices with respect to the base reference frame when all
% the joint variables are zero.
g_00 = [R, p_00; ...
        zeros(1,3), 1];
    
g_01 = [R, p_01; ...
        zeros(1,3), 1];
    
    
g_02 = [R, p_02; ...
        zeros(1,3), 1];
    

g_03 = [R, p_03; ...
        zeros(1,3), 1];

g_12 = inverseTransmat(g_01)*g_02;
g_23 = inverseTransmat(g_02)*g_03;

Mlist = cat(3, g_01, g_12, g_23);
Mi = cat(3, g_01, g_02, g_03);

% Specifying the rotation axis 
% Rotation about the z-axis
% omega = [0;0;1];


% Getting the unit twists for the two revolute joints with respect to the
% base reference frame. The unit twists correspond to the screw axis for
% the revolute joint. 
% q_1 and q_2 are any points of the joint reference frame for the joints 1
% and 2.
q1 = p_00;
q2 = [a1;0;0];
q = [q1, q2];
% q2 = [a1*cos(theta(1));a1*sin(theta(1));0];

twist_1 = GetTwist(omega, q1);twist_2 = GetTwist(omega, q2);
Slist = [twist_1, twist_2];

G1 = blkdiag(m1*eye(3), eye(3));
G2 = blkdiag(m2*eye(3), eye(3));
Glist = cat(3, G1, G2);

% theta_double_dot = [1;0];
% [tau, F_i, V_i, Vdot_i, A_i, Ad_gi] = InverseDynamics(Mi, Mlist, Glist, Slist, theta, theta_dot, theta_double_dot, Ftip, g, q);
% disp(tau);

M = getMassMatrix(Mi, Mlist, Glist, Slist, theta, q);
M
[theta_double_dot_new, Mass_Matrix, C_vector, G_vector, JTFtip] = ForwardDynamics(Mi, Mlist, Glist, Slist, theta, theta_dot, Ftip, tau, g, q);
Mass_Matrix
C_vector
G_vector