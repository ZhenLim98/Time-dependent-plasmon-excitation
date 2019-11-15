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
eF = 1.5*ee;%2*ee;
gamma = 1e12;


alphaF = 0.2;%0.15;
if alphaF > 0.5
    error('Invalid alpha');
end

L = 1; %1e-2
g = pi/L;%0.001*pi/L;
Omega =2*pi/(2.1e-12); % 2*pi*1e12;

Ng = 400;

nqPoints = 1;
qmin = 0.078e6;%0.1e6;
qmax = 0.078e6;%0.1e6;
qVec = linspace(qmin,qmax,nqPoints);

nOmPoints = 1;
omin = 4.4*2*pi*1e12;
omax = 4.4*2*pi*1e12;%2*pi*10.001e12;
omVec = linspace(omin,omax,nOmPoints);

ttVec = zeros(nOmPoints,nqPoints);
rrVec = zeros(nOmPoints,nqPoints);
% Omega = 2*pi*0.200e12;
recVecs = -Ng:Ng;

figure() % %figure()

for jq = 1:length(qVec)
    for jom = 1:length(omVec)
        om = omVec(jom);
        k0 = sqrt(e2)*om/cl;
%         q  = k0*sin(thetaVec(jq));
        q = qVec(jq);
        kz = (q<=k0)*sqrt(k0^2-q^2) + (q>k0)*1i*sqrt(q^2-k0^2);%-(q>k0)*1i*sqrt(q^2-k0^2); % k0*cos(theta);
        if imag(kz)~=0
            warning('evanescent source used')
        end
        check_subwavelength(g,k0);
        sou = zeros(2*Ng+1,1);
        sou(Ng+1) = 1i*2*e2/kz;

         [transmitE,kap1n,kap2n] = trans_calc(q,om,g,Omega,e1,e2,Ng,sou);        
        reflE = 1-transmitE;
        ttVec = e1/e2*kz./abs(kap1n').*abs(transmitE).^2; % n=0 only
        rrVec = (kz./abs(kap1n')).*abs(reflE).^2;
    end
%     hold on
end
% contourf(qVec/1e6,omVec/2/pi/1e12,log(abs(ttVec)),200,'LineColor','none')
    plot((om + recVecs * Omega)/(2*pi*1e12),log(abs(ttVec)),'LineWidth',1.5);
%ylim([0,1])
%xlim([0,10])

% plot(qVec,abs(ttVec(1,:)));

% nlayers = 400;
% figure()
% contourf(thetaVec,omVec/2/pi,abs(ttVec),nlayers,'LineColor', 'none')
% title(['Transmittance, \Omega = ',num2str(Omega)])
% % hold on
% % plot(qVec,Omega/g*qVec,'r')
% colorbar
% ylabel('f')
% xlabel(' \theta   ')
% 
% figure()
% contourf(thetaVec,omVec,abs(rrVec),nlayers,'LineColor', 'none')
% % hold on
% % plot(qVec,Omega/g*qVec,'r')
% colorbar
% title('Reflectance')
% ylabel('$$ \omega $$')
% xlabel('q')
% 
% absVec = 1-ttVec-rrVec;
% 
% figure()
% contourf(thetaVec,omVec,abs(absVec),nlayers,'LineColor', 'none')
% % hold on
% % plot(qVec,omGSPAnVec,'r')
% % hold on
% % plot(qVec,omGSPAnVec_p1,'r')
% % hold on
% % plot(qVec,omGSPAnVec_m1,'r')
% colorbar
% title('Absorbance')
% ylabel('$$ \omega $$')
% xlabel('$$ \theta $$')

% end