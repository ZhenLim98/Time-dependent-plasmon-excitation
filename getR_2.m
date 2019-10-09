
function R = getR_2(k, theta, k2, E_vector, Eincx, n)
%using the BCs from page 171
    R = 0;
    if isreal(k2(n + 1)) == 0
    R = R + (k * cos(theta)/abs(k2(n + 1))) * abs((E_vector(n + 1) - Eincx)/Eincx) * abs((E_vector(n + 1) - Eincx)/Eincx);
    end
    
    for m = 1:1:(n)
        if isreal(k2(m)) == 0
        else
        factor1 = k * ( cos(theta)/ abs(k2(m + n + 1)) ) * abs(E_vector(m + n + 1)/Eincx)* abs(E_vector(m + n + 1)/Eincx);
        factor2 = k * ( cos(theta)/ abs(k2(n + 1 - m)) ) * abs(E_vector(n + 1 - m)/Eincx)* abs(E_vector(n + 1 - m)/Eincx);
        R = R + factor1 + factor2;
        end
    end

end
