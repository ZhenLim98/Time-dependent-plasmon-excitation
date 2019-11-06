function [ out ] = lambda( om_vec,om, Omega, n_vec, Ng)
%conductivity coefficients for a conductivity function that:
% starts at y1 for t = 0
% reaches y2 for t = t1
% goes back to y1 for t = t2
%period length t2
%y axis units is in sigmaDrude
T = 2.1e-12;
y1 = 0.2 ;%* sigmaDrude(om_vec) ;
y2 = 0.8 ;%* sigmaDrude(om_vec) ;
t1 = 0.1e-12;
exp_t = exp(-2* pi * 1i .* t1 .* n_vec ./ T) ;
out = (y1-y2)./(2*pi*1i .* n_vec) .* (1 - 1/(T - t1) * (T - t1 .* exp_t));
out = out + T.*T./(2 * pi .* n_vec .* n_vec .* t1 .* (T - t1)) .* (y1-y2).*(1-exp_t);

%Replaces singularities in the conductivity with the n = 0 fourier coefficient
%change this
nan = isnan(out);
out(nan) =  (y1 + (y2-y1) * T/(2*t1));%(y1(nan) + (y2(nan)-y1(nan)) * T/(2*t1));

%%
%gaussian modulation of the function in frequency space to justify sharp
%cutoff
g_sigma = 500;
mu = 0;
%68% of the content is within +- g sigma here
gaussian =  1/sqrt (2*pi*g_sigma*g_sigma)  * exp(- (n_vec - mu) .* (n_vec-mu)./(2 * g_sigma * g_sigma));
out = out .* gaussian;

end

