%% DYNAMICS OF A 2-R MANIPULATOR USING THE NEWTON-EULER RECURSIVE ALGORITHM.
% By: Aditya Patankar 

clear all;
close all;
clc;

L_1 = 20;
L_2 = L_1;

m1 = 10;
m2 = m1;

% Getting the configuration of the current frame expressed in terms of the
% next frame.

R_0 = [1,0,0;
    0,1,0;
    0,0,1];

p_00 = zeros(3,1);
p_01 = [10;0;0];
p_02 = [30;0;0];

% For the end-effector frame. The reference frame at the tip is not
% oriented as the regular contact frame convention with z-axis pointing
% inward. Its orientation is same as the link and the joint reference
% frames. 
p_03 = [40;0;0];

q_01 = p_01;
q_02 = p_02;
q_03 = p_03;

theta_1 = deg2rad(30);
theta_2 = deg2rad(30);
% For the end effector. 
theta_3 = deg2rad(0);

% Defining the theta vector or the vector of joint angles
theta = [theta_1; theta_2];
% Defining the vector of joint velocities and the vector of joint
% accelerations. 
theta_dot = [4;4];
theta_double_dot = [2;2];

omega = [0;0;1];
omega_3 = zeros(3,1);

g_00 = [R_0, p_00;
    zeros(1,3),1];

g_01 = [R_0, p_01;
    zeros(1,3), 1];

g_02 = [R_0, p_02;
    zeros(1,3), 1];

g_03 = [R_0, p_03;
    zeros(1,3), 1];
 
g_10 = inv(g_01)*g_00;
g_21 = inv(g_02)*g_01;
g_12 = inv(g_01)*g_02;
g_23 = inv(g_02)*g_03;

% Defining the twists and twist_dot along with forces at the tips of the
% end-effectors.
omega_dot = zeros(3,1);
v_dot = [0;-9.8;0];
twist_0 = zeros(6,1);
twist_0_dot = [omega;v_dot];
twist_1 = GetTwist(omega, q_01);
twist_2 = GetTwist(omega, q_02);
F_tip = [10;0;0;0;0;0];

Adjoint_g01 = GetAdjoint(inv(g_01));
Adjoint_g02 = GetAdjoint(inv(g_02));

A_1 = Adjoint_g01*twist_1;
A_2 = Adjoint_g02*twist_2;

exp_twist_theta_1 = GetExponential(omega, theta_1, q_01);
exp_twist_theta_2 = GetExponential(omega, theta_2, q_02);
exp_twist_theta_3 = GetExponential(omega_3, theta_3, q_03);

% Configuration of the reference frame {i} expressed in the reference frame
% {i-1} when we are given the joint variables(theta). Analogous to the body
% form of the product of exponentials formulation.

g_01_theta = g_01*exp_twist_theta_1;
g_12_theta = g_12*exp_twist_theta_2;
g_23_theta = g_23*exp_twist_theta_3;

g_10_theta = inv(g_01_theta);
g_21_theta = inv(g_12_theta);
g_32_theta = inv(g_23_theta);

% This last transformation is required for constructing the L matrix.
g_31_theta = g_32_theta*g_21_theta;

% ad_twist_1 = getSmallAdjoint(twist_1);

% Adjoint of the transformations
Adjoint_g_10_theta = GetAdjoint(g_10_theta);
Adjoint_g_21_theta = GetAdjoint(g_21_theta);
Adjoint_g_32_theta = GetAdjoint(g_32_theta);
Adjoint_g_31_theta = GetAdjoint(g_31_theta);

% Rotational inertia matrices for the links:

I_1 = [0,0,0;
       0, ((1/12)*m1*L_1^2), 0;
       0, 0, ((1/12)*m1*L_1^2)];
   
I_2 = [0,0,0;
       0, ((1/12)*m2*L_2^2), 0;
       0, 0, ((1/12)*m2*L_2^2)];

% Spatial inertia matrix for the links, the inertia due to the joints are
% not considered. 

G_1 = [m1*eye(3,3), zeros(3,3);
      zeros(3,3), I_1];

G_2 = [m2*eye(3,3), zeros(3,3);
      zeros(3,3), I_2];
  
% G matrix

G = [G_1, zeros(6,6);
     zeros(6,6), G_2];

% W matrix 
W_theta = [zeros(6,6), zeros(6,6);
           zeros(6,6), Adjoint_g_21_theta];
       
% A matrix consisting of the transformed twists
A = [A_1, zeros(6,1);
    zeros(6,1), A_2];

% Defining the ad_A_theta_dot matrix:
ad_A1_theta_dot_1 = getSmallAdjoint(A_1*theta_dot(1));
ad_A1_theta_dot_2 = getSmallAdjoint(A_2*theta_dot(2));

ad_A_theta_dot = [ad_A1_theta_dot_1,zeros(6,6);
                  zeros(6,6),ad_A1_theta_dot_2];
              
% Defining some base quantities for the recursive algorithm.
twist_base = [Adjoint_g_10_theta*zeros(6,1); zeros(6,1)];
twist_dot_base = [Adjoint_g_10_theta*twist_0_dot; zeros(6,1)];
F_tip_transformed = [zeros(6,1);Adjoint_g_32_theta'*F_tip];

% Defining the L matrix for 2R manipulator:
L_theta = [eye(6,6), zeros(6,6);
           Adjoint_g_21_theta, eye(6,6)];


%% Computing the quantities using the recursive algorithm

% Computing the twists
twist = L_theta*(A*theta_dot + twist_base);

% Computing the accelerations
twist_dot = L_theta*(A*theta_double_dot + ad_A_theta_dot*W_theta*twist + ad_A_theta_dot*twist_base + twist_dot_base);

% Computing the ad_twist matrix required for computing the forces.
 computed_twist_1 = twist(1:6);
 computed_twist_2 = twist(7:12);
 
 ad_twist_1 = getSmallAdjoint(computed_twist_1);
 ad_twist_2 = getSmallAdjoint(computed_twist_2);
 
 ad_twist = [ad_twist_1, zeros(6,6);
             zeros(6,6), ad_twist_2];
         
% Computing the forces generated:
computed_F = L_theta'*(G*twist_dot - ad_twist'*G*twist + F_tip_transformed);

% Computing the generated torques:
tau = A'*computed_F;

% Computing the mass matrix, Coriolis matrix and the gravitational effect
% vector:

M = A'*L_theta'*G*L_theta*A;

C = -A'*L_theta'*(G*L_theta*ad_A_theta_dot*W_theta + ad_twist'*G)*L_theta*A*theta_dot;

g = A'*L_theta'*G*L_theta*twist_dot_base;
 
 
 