function [tau, F_i, V_i, Vdot_i] = InverseDynamics(Mlist, Glist, Slist, theta, theta_dot, theta_double_dot, Ftip)

    n = size(theta,1);

    % Initializing the multidimensional arrays
    g = [0;0;9.8];
    % Initializing Mi to a 4*4 identity matrix is incorrect and this needs 
    % be rectified.
    Mi = eye(4);
    A_i = zeros(6,n);

    % Since the last matrix contains the transformation between the last link
    % reference frame and the end-effector reference frame.
    Ad_gi = zeros(6,6,n+1);

    % Computed twists and accelerations.
    V_i = zeros(6, n+1);

    % Initializing the acceleration for the base case.;
    Vdot_i = zeros(6, n+1);
    Vdot_i(1:3,1) = -g;

    F_i = zeros(6, n+1);
    F_i(:,n+1) = Ftip;
    tau = zeros(n,1);

    Ad_gi(:,:,n+1) = GetAdjoint(inverseTransmat(Mlist(:,:,n+1)));
    % Forward iterations of the Newton-Euler inverse dynamics algorithm.
    for i = 1:n
        A_i(:,i) = GetAdjoint(inverseTransmat(Mi))*Slist(:,i);
        Ad_gi(:,:,i) = GetAdjoint(se3toSE3(Slist(:,i), -theta(i))*inverseTransmat(Mlist(:,:,i)));
        V_i(:,i+1) = A_i(:,i)*theta_dot(i) + Ad_gi(:,:,i)*V_i(:,i);
        Vdot_i(:,i+1) = A_i(:,i)*theta_double_dot(i) + Ad_gi(:,:,i)*Vdot_i(:,i) + getSmallAdjoint(V_i(:,i+1))*A_i(:,i)*theta_dot(i);
    end

    % Backward iterations of the Newton-Euler inverse dynamics algorithm.

    for i = n:-1:1
       F_i(:,i) = Ad_gi(:,:,i+1)*F_i(:,i+1) +  Glist(:,:,i)*Vdot_i(:,i+1) - (getSmallAdjoint(V_i(:,i+1)))'*Glist(:,:,i)*Vdot_i(:,i+1);
       tau(i) = F_i(:,i)'*A_i(:,i);
    end


end
