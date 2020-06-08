function ad_twist = getSmallAdjoint(twist)

v = twist(1:3);
omega = twist(4:6);

v_hat = skewSymmetric(v);
omega_hat = skewSymmetric(omega);

ad_twist = [omega_hat, v_hat;
            zeros(3,3), omega_hat];

end 