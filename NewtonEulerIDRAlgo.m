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

g = [0;0;9.8];

% Computing the torques and forces required at each joints using the
% Newton-Euler Recursive Algorithm for inverse dynamics.
[tau, F_i, V_i, Vdot_i, A_i, Ad_gi] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, theta_double_dot, Ftip, g);



