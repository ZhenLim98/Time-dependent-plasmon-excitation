clear all
global hbar ee cl eF gamma e0

% for Omega = 2*pi*linspace(0,0.4,5)*1e12

hbar = 1.05e-34;
ee = 1.6e-19;
cl = 3e8;
e0 = 8.85e-12;
fsConst = 1/137;
e1 = 1;
e2 = 1;
gamma = 1e11;

x = linspace(1,10,10);
y = linspace(1,10,10);

[X,Y] = meshgrid(x,y);

% Z = X.*Y;
% figure()
% contourf(X,Y,Z)

nom = 200;
nk = 200;

kVec = linspace(1e4,2e6,nk);
omVec = 2*pi*linspace(1e11,10e12,nom);
[K,OM] = meshgrid(kVec,omVec);

figure()
for eF = [2]*ee
    DrW = ee^2/pi/hbar^2*eF;
    disp_vis = e1./sqrt(K.^2-e1*OM.^2/cl^2) ...
        + e2./sqrt(K.^2-e2*OM.^2/cl^2) ...
        - DrW/e0./(OM.^2);
    contourf(K/1e6,OM/2/pi/1e12,log(abs(real((1./disp_vis)))),200,'LineColor','none')
    hold on
%     plot(kVec/1e6,sqrt(DrW*kVec/(e1+e2)/e0)/2/pi/1e12,'LineWidth',2.5)
end

% caxis([1e-5,1e-2])
% xlabel('rad/\mum')
% ylabel('f [THz]')
ylim([0.1,10])
xlim([0.01,0.5])
% title('Full')