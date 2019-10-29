function [ out ] = sigma2( om_vec,om, Omega, n_vec, Ng)
%conductivity coefficients for a conductivity function that:
% starts at y1 for t = 0
% reaches y2 for t = t1
% goes back to y1 for t = t2
%period length t2
%y axis units is in sigmaDrude
sigma_drude = sigmaDrude(om_vec);
T = 2 * pi/Omega;
y1 = 0.2;
y2 = 0.8;
t1 = 0.1;
t2 = 2.1;
grad1 = 6;
grad2 = -0.6/2;

out = sigma_drude./(2.*T) .* ((y1 .* (exp( -1i .* n_vec )) - 1)./(-1i .* n_vec .* Omega) + ...
                        grad1 .* sigma_drude./(2 .* T) .* (1./(n_vec .* n_vec .* Omega .* Omega))...
                        .* (exp(-1i .* n_vec .* Omega .* t1) - 1) + ...
                        grad1 .* sigma_drude./(2 .* T).*(t1.* exp(- n_vec .* Omega .* t1))/(-1i .* n_vec .* Omega));
%first section of fourier expansion

out = out + sigma_drude./(2*T) * ((y2)./(-1i .* n_vec .* Omega) ...
                        .* (exp(-1i .* n_vec .* t2) - exp(-1i .* n_vec .* Omega .* t1))...
                         + grad2 .*(sigma_drude./(2 * T))./(-1i .* n_vec .* Omega) ...
                     .* (t2 .* exp(-1i .* n_vec .* Omega .* t2) - t1*(exp(-1i .* n_vec .* Omega .* t1)))...
                     + (1./(n_vec .* n_vec .* Omega .* Omega).*(exp(-1i .* n_vec .* Omega .* t2) - (exp(-1i .* n_vec .* Omega .* t1)))));
%second section of the fourier expansion
%for the n = 0 fourier expansion coefficient.
sigma_0 = sigmaDrude(om);
out(Ng + 1) = sigma_0/(2 .* T) .* (y1 .* t1 + grad1 .* t1 .* t1 ./ 2)  ...
                        + sigma_0/(2 .* T) .* ((t2 - t1) .* y2 + grad2./2 .* (t2 .* t2 - t1 .* t1) ); 
                    
%CONVERTING CONDUCTIVITY INTO DRUDE WEIGHT:
% out = out * -1i * omega/c
end

