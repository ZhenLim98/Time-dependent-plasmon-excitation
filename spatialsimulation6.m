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

c = 3e8;
eps_0 = 8.85e-12;
h = 6.63e-34;
e_e = 1.6 * (10e-19);
sigma_0 = pi * e_e ^ 2 /(2 * h) ;
G = 2 * pi / R;
nSteps = round((f_top - f_bot)/spacing);

Tdata = linspace(0.0, 0.0, nSteps);
f = linspace(0.0, 0.0, nSteps);
Rdata = linspace(0.0, 0.0, nSteps);

% x = sym('x');

for frequency = f_bot : spacing : f_top
    omega = 2 * pi * frequency;
    %calculating R, T over a given range of frequencies
    j1 = round((frequency - f_bot)/spacing) + 1;
    k =  sqrt(e_2) * omega/c;
    q = k * sin(theta);
  
    [k1, k2] = getk_2(n, G, q, e_1, e_2, omega, c);
    
    drude_conductivity =  (4 * E_f) * (sigma_0/(pi * (h / (2 * pi))))  * (1/ ( (gamma) - 1i * omega));
    M = getM_2(n, drude_conductivity, d_g, R, omega, eps_0, e_1, e_2, k1, k2);
    
    Einit = zeros(1 , 2 * n + 1)';
    Einit(n + 1) = 1i * 2 * e_2 * Eincx /(k * cos(theta)) ;
    E_vector = M \ (Einit);
    
    Tdata(j1) = getT_2(n, k, theta, k1, Eincx, E_vector, e_1, e_2);
    Rdata(j1) = getR_2(k, theta, k2, E_vector, Eincx, n);
    
    f(j1) = frequency;
end    
end
