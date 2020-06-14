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
[tau, F_i, V_i, Vdot_i, A_i, Ad_gi] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, theta_double_dot, Ftip,g);

%% Computing the Dynamics equation for 2R manipulator in closed form:

m = size(A_i, 2);
n = size(V_i,1);

% Initializing the empty stacked vectors and block diagonal matrices.
V = zeros(6*m,1);
F = zeros(6*m,1);
V_base = zeros(6*m,1);
Vdot_base = zeros(6*m,1);
F_stacked = zeros(6*m,1);
ad_Vi = zeros(6,6,m);
temp_Ad_gi = zeros(6,6,m);
ad_A_theta_dot = zeros(6,6,m);


% Concatenating the V_i's to construct the vertically stacked vector 'V'
for i = 1:m
    if i == 1
        V(1:n,:) = V_i(:,i+1);
    elseif i > 1
        V(((i-1)*(n+1)):(i*n),:) = V_i(:,i+1);
    end
end

%  Concatenating the F_i's to construct the vertically stacked vector 'F'
for i = 1:m
    if i == 1
        F(1:n,:) = F_i(:,i);
    elseif i > 1
        F(((i-1)*(n+1)):(i*n),:) = F_i(:,i);
    end
end

% Constructing the block diagonal matrix A for storing the transformed unit
% twists A_i's.
temp_mat = mat2cell(A_i, n, [1,1]);
A = blkdiag(temp_mat{:});

% Constructing the block diagonal matrix G for storing the spatial inertia
% matrices G_i's.
temp_mat = mat2cell(Glist, n, n, [1,1]);
G = blkdiag(temp_mat{:});

% Adding elements to the ad_Vi matrix using a for loop.
for i = 1:m
    ad_Vi(:,:,i) = getSmallAdjoint(V_i(:,i+1));
end
% Constructing the block diagonal matrix for storing all the small adjoint
% of the computed twists.
temp_mat = mat2cell(ad_Vi, 6, 6, [1,1]);
ad_Vi = blkdiag(temp_mat{:});


% Constructing the block diagonal matrix for small adjoint of the product
% of A_i and theta dot.
for i = 1:m
   ad_A_theta_dot(:,:,i) = getSmallAdjoint(A_i(:,i)*theta_dot(i)); 
end
temp_mat = mat2cell(ad_A_theta_dot, n, n, [1,1]);
ad_A_theta_dot = blkdiag(temp_mat{:});


% Constructing the W matrix using the Adjoint of the transformation
% matrices.
temp_Ad_gi(:,:,2) = Ad_gi(:,:,2);
temp_mat = mat2cell(temp_Ad_gi, n, n, [1,1]);
W = blkdiag(temp_mat{:});


% Constructing the L matrix for a 2R manipulator.
L = cat(3,eye(6,6),eye(6,6));
temp_mat = mat2cell(L,6,6,[1,1]);
L = blkdiag(temp_mat{:});
L(7:12, 1:6) = Ad_gi(:,:,2);


% Computing the stacked vectors denoting the initial conditions such as
% V_base, Vdot_base and F_stacked.
V_base(1:6,:) = Ad_gi(:,:,1)*V_i(:,1);
Vdot_base(1:6,:) = Ad_gi(:,:,1)*Vdot_i(:,1);
F_stacked(7:12,:)= Ad_gi(:,:,3)'*F_i(:,3);

%% Computing the terms of the Dynamics Equation for a 2R Manipulator:

Mass_Matrix = A'*L'*G*L*A;
C_Vector = -A'*L'*(G*L*ad_A_theta_dot*W + ad_Vi'*G)*L*A*theta_dot;
G_Vector = A'*L'*G*L*Vdot_base;
