close all;
clear all;
clc;

% Reading the required data from the text file.
data = getData('2R.txt');
[x,~] = size(data{1});
C = {};

% In the for loop we assign a value to a variable which is then used for
% our processing purposes.
% The variable names are stored in the first column of the cell array 'C'
% whereas the values to be assigned are stored in the third column.

for i=1:x
% Each element of the cell 'data' is split and stores in another cell
% array of size m*n.
% Each row of 'C' contains the split string.
    C{i} = strsplit(data{1}{i}, ' ');
    assignin('base', C{i}{1}, str2num(C{i}{3}));
end

% Defining the transformation matrices for the individual link reference
% frames with respect to the base reference frame. 

g_00 = [R_0, p_00;
    zeros(1,3),1];

g_01 = [R_0, p_01;
    zeros(1,3), 1];

g_02 = [R_0, p_02;
    zeros(1,3), 1];

g_03 = [R_0, p_03;
    zeros(1,3), 1];

g_12 = inverseTransmat(g_01)*g_02;
g_23 = inverseTransmat(g_02)*g_03;

% Concatenating all the transformation matrices in a single
% multidimensional array.
Mi = cat(3, g_01, g_02, g_03);
Mlist = cat(3, g_01, g_12, g_23);

% Concatenating the spatial inertial matrix into a single multidimensional
% array.
G1 = diag(g1);
G2 = diag(g2);

Glist = cat(3, G1, G2);

% Computing the twists and concatenating them into a single
% multidimensional array
q_1 = p_01;
q_2 = p_02;

twist_1 = GetTwist(omega, q_1);
twist_2 = GetTwist(omega, q_2);

twist_list = [twist_1, twist_2];

% This is the input to our Forward Dynamics function. 
tau = [0.5; 0.6];
g = [0;9.8;0];

% theta_double_dot_new = FowardDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, Ftip, tau, g);

%% Euler Integration code.
dt = 0.1;
Res = 8;
n = size(tau,1);
% Random values of the joint torques/efforts required at each time step.
taumat = [[3.63; -6.58], [3.74; -5.5], ...
        [4.31; -5.19], [5.18; -4.31], ...
        [5.85; -2.59], [5.78; -1.7], ...
        [4.99; -1.19], [4.08; 0.07], ...
        [3.56; 0.97], [3.49; 1.23]];

% Random values of the end-effector forces.
Ftip = ones(6,size(taumat,2));

% Initializing the theta 
theta_new_mat = zeros(n,size(taumat,2));
theta_new_mat(:,1) = theta;
theta_dot_new_mat = zeros(n, size(taumat,2));
theta_dot_new_mat(:,1) = theta_dot;
theta_double_dot_mat = zeros(n, size(taumat,2));

for i = 1:size(taumat,2)-1
    for j = 1:Res
        theta_double_dot_new = FowardDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, Ftip(:,i), taumat(:,i), g);
        theta_double_dot_mat(:,i) = theta_double_dot_new;
        alpha= dt/Res;
        [theta_next, theta_dot_next] = EulerIntegration(theta_new_mat(:,i), theta_dot_new_mat(:,i), theta_double_dot_new, alpha);
    end
    theta_new_mat(:,i+1) = theta_next;
    theta_dot_new_mat(:, i+1) = theta_dot_next;
end

%% Plotting the joint space variables 

Tf = size(taumat, 1);
time=0: (Tf / size(theta_new_mat, 2)): (Tf - (Tf / size(theta_new_mat, 2)));
plot(time,theta_new_mat(1, :),'b')
hold on
plot(time,theta_new_mat(2, :), 'g')
plot(time,theta_dot_new_mat(1, :), 'c')
plot(time,theta_dot_new_mat(2, :), 'm')
title('Plot of Joint Angles and Joint Velocities')
xlabel('Time')
ylabel('Joint Angles/Velocities')
legend('Theta1', 'Theta2', 'Theta3', 'DTheta1', 'DTheta2', 'DTheta3')
