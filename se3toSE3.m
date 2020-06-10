function exp_twist_hat_theta = se3toSE3(twist, theta)


% In this function, the twist, which is one of the input arguments, contains the linear velocity of translation
% and the angular velocity of rotation. 

    I = eye(3);
    
%     Omega corresponds to the unit rotation axis, and theta corresponds to
%     the angle of rotation about that axis.
    v = twist(1:3);
    omega = twist(4:6);
    
    exp_omega_hat_theta = so3toSO3(omega, theta);
    
%     When angular velocity is not equal to zero.
    
    exp_twist_hat_theta = [exp_omega_hat_theta, (I - exp_omega_hat_theta)*cross(omega, v) + omega*omega'*v*theta;
                           zeros(1,3),           1];
    
    
end