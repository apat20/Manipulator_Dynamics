function theta_double_dot_new = FowardDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, Ftip, tau, g)

    % Computing the mass matrix for use in Euler integration of manipulator
    % equations of motion.
    n = size(theta,1);
    g_temp = [0;0;0];
    Ftip_temp = [0;0;0;0;0;0];
    theta_dot_temp = zeros(size(theta,1),1);
    Mass_Matrix = zeros(n);

    for i = 1:n
       theta_double_dot_temp =  zeros(size(theta,1),1);
       theta_double_dot_temp(i) = 1;
       [tau, ~, ~, ~, ~, ~] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot_temp, theta_double_dot_temp, Ftip_temp, g_temp);
       Mass_Matrix(:,i) = tau;
    end

    % Computing the Coriolis effect vector for use in Euler integration of
    % manipulator equations of motion.
    theta_double_dot_temp = zeros(size(theta,1),1);
    [C_vector, ~, ~, ~, ~, ~] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot, theta_double_dot_temp, Ftip_temp, g_temp);


    % Computing the gravity vector for use in the Forward dynamics equation.
    [G_vector, ~, ~, ~, ~, ~] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot_temp, theta_double_dot_temp, Ftip_temp, g);


    % Computing the joint forces and torques required to create the force at
    % the end-effector Ftip. 
    [JTFtip, ~, ~, ~, ~, ~] = InverseDynamics(Mi, Mlist, Glist, twist_list, theta, theta_dot_temp, theta_double_dot_temp, Ftip, g_temp);


    % Using the Forward Dynamics equation to compute the values of theta_double
    %_dot for the computed values of the Mass_matrix, Coriolis vector, Gravity vector
    % and the Torques required for the end effector forces given a vector of
    % torques, joint angles and joint velocities using the Backward slash operator. 
    theta_double_dot_new = Mass_Matrix\(tau - C_vector - G_vector - JTFtip);
    
end