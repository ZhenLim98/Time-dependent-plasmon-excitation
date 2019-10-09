
    function T = getT_2(n, k, theta, k1, Eincx, E_vector, e_1, e_2)
        %using the BCs from page 171
        T = 0;
        
        for j = 1:1:(2 * n + 1)
            %if k1 is real then the wave is decaying
            %R and T only include propagating modes
            if isreal(k1(j)) == 0 
            else
            factor1 = (e_1 / e_2 ) * (k * cos(theta) / abs(k1(j))) * abs(E_vector(j)/Eincx) * abs(E_vector(j)/Eincx);
            T = T + factor1;
            end
        end
    end
