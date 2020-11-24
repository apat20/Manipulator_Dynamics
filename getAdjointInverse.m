function Adjoint_inverse = getAdjointInverse(g)
    R = g(1:3,1:3);
    p = g(1:3,4);
    
    p_hat = skewSymmetric(p);
    
    Adjoint_inverse = [R', -R'*p_hat; ...
                        zeros(3), R'];
end