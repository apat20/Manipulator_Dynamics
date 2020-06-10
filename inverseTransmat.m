% This function is used to compute the inverse of a transformation matrix 


function inv_g = inverseTransmat(g)
    R = g(1:3,1:3);
    p = g(1:3,4);
    
    inv_g = [R', -R'*p;
            zeros(1,3), 1];
        
end