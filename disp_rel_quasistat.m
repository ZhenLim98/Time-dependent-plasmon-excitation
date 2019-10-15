clear all
global hbar ee cl eF e0

% for Omega = 2*pi*linspace(0,0.4,5)*1e12

hbar = 1.05e-34;
ee = 1.6e-19;
cl = 3e8;
e0 = 8.85e-12;
fsConst = 1/137;
e1 = 1;
e2 = 1;
gamma = 1e11;

nk = 100;
kVec = linspace(5e4,1e6,nk);


figure()
for eF = [0.2,0.6,1.0,1.4,1.8]*ee
    DrW = ee^2/pi/hbar^2*eF;
    omVec = sqrt(DrW*kVec/(e1+e2)/e0);
    plot(kVec/1e6,omVec/2/pi/1e12,'LineWidth',2.5)
    hold on
end

plot(kVec/1e6,cl*kVec/2/pi/1e12,'--k','LineWidth',2.5)
xlabel('rad/\mum')
ylabel('f [THz]')
ylim([0.1,15])
title('Quasistatic')




