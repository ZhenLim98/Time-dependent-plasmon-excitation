
    function M = getM_2(n, drude_conductivity, d_g, R, omega, eps_0, e_1, e_2, k1, k2 )
        % in the case of a physical grating
        %Drude conductivity used (~THz regime)
        %finding the fourier coefficients of the conductivity as a vector
        
        %constructing conductivity coefficients
        l_vector = (-n : n);
        sigma_n = (drude_conductivity ./(pi .* l_vector)) .* (sin(l_vector .* pi * d_g/R) .* exp(1i .* l_vector * d_g * pi/R));
        sigma_n(n + 1) = drude_conductivity * d_g/R;

        %constructing diagonals of matrix M (conductivity sum)
        diagmain = repelem(sigma_n(n + 1), 2 * n + 1);
        conductivity_matrix = zeros(2 * n + 1, 2 * n + 1);
        
        for diag_index = 1:1:(n)
            a = repelem(sigma_n(n + diag_index + 1), 2 * n - diag_index + 1);
            b = repelem(sigma_n(n + 1 - diag_index), 2 * n - diag_index + 1);
            conductivity_matrix = conductivity_matrix + diag(a, +diag_index) + diag(b, -diag_index);
        end
        
        conductivity_matrix = conductivity_matrix + diag(diagmain);
        M = ( 1i/(omega * eps_0) ) .* conductivity_matrix;
        
        %constructing additional diagonal term of M (first term on the left)
        
        diagonalterm = (e_1./(k1) + e_2./(k2));
        diagonal = diag(diagonalterm);
        M = M + diagonal;
       
    end