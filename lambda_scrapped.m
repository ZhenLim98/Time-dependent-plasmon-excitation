function [ out ] = lambda(  n_vec)
%conductivity coefficients for a conductivity function that:
% starts at y1 for t = 0
% reaches y2 for t = t1
% goes back to y1 for t = t2
%period length t2
%y axis units is in sigmaDrude
T = 2.1e-12;
Drude_w = 76e9; 
y1 = 0.2 * Drude_w ;
y2 = 0.8*  Drude_w;
t1 = 0.1e-12;
exp_t = exp(-2* pi * 1i .* t1 .* n_vec ./ T) ;
out = (y1-y2)./(2*pi*1i .* n_vec) .* (1 - 1/(T - t1) * (T - t1 .* exp_t));
out = out + T.*T./(2 * pi .* n_vec .* n_vec .* t1 .* (T - t1)) .* (y1-y2).*(1-exp_t);

%Replaces singularities in the conductivity with the n = 0 fourier coefficient
%change this
zero = find(~n_vec);
out(zero) =  (y1 + (y2-y1) * T/(2*t1));

%%
%gaussian modulation of the function in frequency space to justify sharp
%cutoff
g_sigma = 4;
mu = 0;
%68% of the content is within +- g sigma here
gaussian =  1/sqrt (2*pi*g_sigma*g_sigma)  * exp(- abs(n_vec - mu) .* abs(n_vec-mu)./(2 * g_sigma * g_sigma));
out = out .* gaussian;
% %% Finding and scaling IFT
% L=length(out);                      % series length
% Nsamp = 201; 
% Fs=Nsamp/(T*1e12);  % 3 points (y1, y2 then y1) per T
% t = (0: 10/Fs: 1 - 10/Fs) * 2.1e-12;
% ift = abs(ifft(out)) * L;
% % ift = repmat(ift, 1,Fs/L);
% plot(t, ift)

end

