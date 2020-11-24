%This function is used to implement the exponential of twist formula 
%The exponential of twists formula is used to calculate the rigid body
%transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%INPUT ARGUMENTS%%

% 'Omega' is a 3*1 vector which which encodes the information regarding the axis 
%rotation.

% 'q' is a 3*1 vector which contains the position of a point on the axis of
%rotation. The choice of this point is arbitary and depends on the user.

% 'theta' is the angle of rotation in degrees. The input is in radians
% itself.


%%
function exp_twist_theta = GetExponential(twist, theta, q)
    
    %The identity matrix 'I'
    I = eye(3);
    
%     Omega corresponds to the unit rotation axis, and theta corresponds to
%     the angle of rotation about that axis.
    v = twist(1:3);
    omega = twist(4:6);
    
    %Skew symmetric form of the omega vector
    omega_hat = skewSymmetric(omega);
    
   %Using the Rodriguez formula to calculate the exponential of omega
   %formula to get the exponential coordinates for rotation.
   
   % Write a function for this and include it in here.!!!!!
   exp_omega_hat_theta = I + omega_hat*sin(theta) + omega_hat^2*(1-cos(theta));
   
   %Calculating the transformation using the exponential of twist formula
   exp_twist_theta = [exp_omega_hat_theta, (I - exp_omega_hat_theta)*q;
                      zeros(1,3),           1];
% % The expression below is derived from Modern Robotics by Kevin Lynch.            
%    exp_twist_theta = [exp_omega_hat_theta, (I*theta + (1-cos(theta))*omega_hat + (theta-sin(theta))*omega_hat^2)*v;
%                       zeros(1,3),           1];               
   
end