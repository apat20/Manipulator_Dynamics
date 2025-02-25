% This function is used to compute the joint angles and the joint rates
% using the Euler step implemented in the Euler first order iteration. We
% use to compute the joint variables(angles and rates) to be used for
% computing the Forward Dynamics of an n-DOF manipulator.


function [theta_next, theta_dot_next] = EulerIntegration(theta_current, theta_dot_current, theta_double_dot, alpha)
    theta_next = theta_current + alpha*theta_dot_current;
    theta_dot_next = theta_dot_current + alpha*theta_double_dot;
end