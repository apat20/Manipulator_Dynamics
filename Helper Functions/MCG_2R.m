function [M, C, G] = MCG_2R(m1, a1, lc1, I1, m2, a2, lc2, I2, g, theta1, theta2, thetad1, thetad2)

M = [...
 m2*a1^2 + 2*m2*cos(theta2)*a1*lc2 + m1*lc1^2 + m2*lc2^2 + I1 + I2, m2*lc2^2 + a1*m2*cos(theta2)*lc2 + I2
                                I2 + lc2*m2*(lc2 + a1*cos(theta2)),                         m2*lc2^2 + I2];

C = [...
 -a1*lc2*m2*thetad2*sin(theta2)*(2*thetad1 + thetad2)
                      a1*lc2*m2*thetad1^2*sin(theta2)];

G = [...
 g*lc2*m2*cos(theta1 + theta2) + a1*g*m2*cos(theta1) + g*lc1*m1*cos(theta1)
                                              g*lc2*m2*cos(theta1 + theta2)];
end