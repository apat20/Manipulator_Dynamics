function exp_omega_hat_theta = so3toSO3(omega, theta)
    %The identity matrix 'I'
    I = eye(3);
    
    %Skew symmetric form of the omega vector
    omega_hat = skewSymmetric(omega);
    
   %Using the Rodriguez formula to calculate the exponential of omega
   %formula to get the exponential coordinates for rotation.
   exp_omega_hat_theta = I + omega_hat*sin(theta) + omega_hat^2*(1-cos(theta));
end