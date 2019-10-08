function [Rdata, Tdata, f] = spatialsimulation6(e_1, e_2, Eincx, R, theta, E_f, gamma, f_bot, f_top, spacing, n, d_g)
%For plane waves incident on a sheet of graphene, calculates the
%Reflectance and Transmittance over a range of frequencies
% e_1           electric permittivity of the SECOND material,
% e_2           electric permittivity of the 1st material
% theta         angle to the vertical
% solving matrix equation M * E_vector = Einit
% Drude conductivity of graphene used - only valid in THz regime
%using Standard Units except energy = eV

%Example command: [Rdata, Tdata, f] = spatialsimulation6(4, 3 ,1, 8 * 10^(-6), 0 , 0.45 * 1.6 * 10^(-19), 3.7 * 10^(-3)  * 1.6 * 10^(-19)/(6.63*10^(-34)/(2*pi)) , 10^(12),  10 * 10^(12) , 10^(11), 100, 4 * 10^(-6))

c = 3 * 10 ^ 8;
eps_0 = 8.85 * 10^(-12);
h = 6.63 * 10^(-34);
e_e = 1.6 * 10^(-19);
sigma_0 = pi * e_e ^ (2)/(2 * h) ;
G = 2 * pi / R;
nSteps = round((f_top - f_bot)/spacing);

Tdata = linspace(0.0, 0.0, nSteps);
f = linspace(0.0, 0.0, nSteps);
Rdata = linspace(0.0, 0.0, nSteps);

x = sym('x');
M = zeros(2 * n + 1 , 2 * n + 1 );
k1 = zeros(1, 2 * n + 1);
k2 = zeros(1, 2 * n + 1);

for frequency = f_bot : spacing : f_top
    omega = 2 * pi * frequency;
    %calculating R, T over a given range of frequencies
    j1 = round((frequency - f_bot)/spacing) + 1;
    k =  sqrt(e_2) * omega/c;
    q = k * sin(theta);
    getk();
    drude_conductivity =  (4 * E_f) * (sigma_0/(pi * (h / (2 * pi))))  * (1/ ( (gamma) - 1i * omega));
    getM();
    Einit = zeros(1 , 2 * n + 1)';
    Einit(n + 1) = 1i * 2 * e_2 * Eincx /(k * cos(theta)) ;
    E_vector = M \ (Einit);
    Tdata(j1) = getT();
    Rdata(j1) = getR();
    f(j1) = frequency;
end

%% Nested functions below]=
    function getk()
        
        G_vector = linspace(-n, n, 2 * n + 1) * G;
        k1 = sqrt( (q +  G_vector) .^ 2 - e_1 * omega * omega / (c * c) );
        k2 = sqrt( (q +  G_vector) .^ 2 - e_2 * omega * omega / (c * c) );
        
    end

    function getM()
        % in the case of a physical grating
        %Drude conductivity used (~THz regime)
        %finding the fourier coefficients of the conductivity as a vector
        
        %constructing conductivity coefficients
        l_vector = linspace(- n , n,  2 * n + 1);
        sigma_n = (drude_conductivity ./(pi .* l_vector)) .* (sin(l_vector .* pi * d_g/R) .* exp(1i .* l_vector * d_g * pi/R));
        sigma_n(n + 1) = drude_conductivity * d_g/R;
%          sigma_n = flip(sigma_n);
        %constructing diagonals of matrix M (conductivity sum)
        diagmain = repelem(sigma_n(n + 1), 2 * n + 1);
        conductivity_matrix = zeros(2 * n + 1, 2 * n + 1);
        for diag_index = 1:1:(n)
            a = repelem(sigma_n(n + diag_index + 1), 2 * n - diag_index + 1);
            b = repelem(sigma_n(n + 1 - diag_index), 2 * n - diag_index + 1);
            conductivity_matrix = conductivity_matrix + diag(a, diag_index) + diag(b, -diag_index);
        end
        
        conductivity_matrix = conductivity_matrix + diag(diagmain);
        M = ( 1i/(omega * eps_0) ) .* conductivity_matrix;
        
        %constructing additional diagonal term of M (first term on the left)
        
        diagonalterm = (e_1./(k1) + e_2./(k2));
        diagonal = diag(diagonalterm);
        M = M + diagonal;
       
    end

    function T = getT()
        %using the BCs from page 171
        
        T = 0;
        for j = 1:1:(2 * n + 1)
            factor1 = (e_1 / e_2 ) * (k * cos(theta) / abs(k1(j))) * abs(E_vector(j)/Eincx) * abs(E_vector(j)/Eincx);
            T = T + factor1;
        end
    end

    function R = getR()
        %using the BCs from page 171
        R = 0;
        R = R + (k * cos(theta)/abs(k2(n + 1))) * abs((E_vector(n + 1) - Eincx)/Eincx) * abs((E_vector(n + 1) - Eincx)/Eincx);
        for m = 1:1:(n)
            factor1 = k * ( cos(theta)/ abs(k2(m + n + 1)) ) * abs(E_vector(m + n + 1)/Eincx)* abs(E_vector(m + n + 1)/Eincx);
            factor2 = k * ( cos(theta)/ abs(k2(n + 1 - m)) ) * abs(E_vector(n + 1 - m)/Eincx)* abs(E_vector(n + 1 - m)/Eincx);
            R = R + factor1 + factor2;
        end
            
    end


end
