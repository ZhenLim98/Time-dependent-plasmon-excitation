function [ out ] = lambda_i(  n_vec,length_vec, j, Ng, Omega)
%%
%Drude weight coefficients for a Drude Weight function that:
% starts at y1 for t = 0
% reaches y2 for t = t1
% goes back to y1 for t = t2
    %period length t2
%This code section could be wrong!

Drude_w = 76.11e9   ;     
y1 = (0.2 )* Drude_w  ;
y2 = (0.8 )*  Drude_w ;
T = 2*pi/Omega; %T is the period of the pump, same as the period of the modulation
t1 = 0.05*T; %rising time

exp_t = exp(-2* pi * 1i .* t1 .* n_vec ./ T) ;
out = (y1-y2)./(2*pi*1i .* n_vec) .* (1 - 1/(T - t1) * (T - t1 .* exp_t));
out = out + T.*T./(2 * pi .* n_vec .* n_vec .* t1 .* (T - t1)) .* (y1-y2).*(1-exp_t);

%Replaces singularities in the conductivity with the n = 0 fourier coefficient
%change this
zero_pos = find(~n_vec);
out(zero_pos) =  (y1 + (y2-y1) * T/(2*t1));

%multiply with gaussian profile in frequency space to justify cutoff
w_falloff = 200;
mu = 0;
%68% of the content is within +- g sigma here
gaussian = exp(- abs(n_vec - mu)./( w_falloff));
out = out .* gaussian;

%%
%repeating the modulation
out = out(Ng + 1 + j);
out = repmat(out, 1, numel(length_vec));
%%

end


% %% Finding and scaling IFT
% L=length(out);                      % series length
% Nsamp = 201; 
% Fs=Nsamp/(T*1e12);  % 3 points (y1, y2 then y1) per T
% t = (0: 10/Fs: 1 - 10/Fs) * 2.1e-12;
% ift = abs(ifft(out)) * L;
% % ift = repmat(ift, 1,Fs/L);
% plot(t, ift)
