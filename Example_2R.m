clc
clearvars

m1 = 1; % Kg
a1 = 0.4; % m
lc1 = a1/2; % m
I1 = 1;
% I1 = 1/12*m1*a1^2;  % Kg.m^2 (Moment of inertia relative to the CoM of the link and about the Z0-axis, refer to the figure)
                    % This is not a matrix here.

m2 = 1; % Kg
a2 = 0.4; % m
lc2 = a2/2; % m
I2 = 1;
% I2 = 1/12*m2*a2^2; % Kg.m^2 (Moment of inertia relative to the CoM of the link and about the Z0-axis, refer to the figure)
                   % This is not a matrix here.

theta1 = 0; % rad
theta2 = pi/2; % rad
thetad1 = 0.2; % rad/s
thetad2 = 0.2; % rad/s

g = 9.8; % m/s^2 (This is just the "magnitude" of gravitational acceleration, for the direction, refere to the figure)

[M, C, G] = MCG_2R(m1, a1, lc1, I1, m2, a2, lc2, I2, g, theta1, theta2, thetad1, thetad2);

disp("M =")
disp(M)

disp("C =")
disp(C)

disp("G =")
disp(G)
