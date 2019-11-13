function [ out,kap1n,kap2n ] = mMat( q,om,g,Omega,e1,e2,Ng,alphaF )
global e0

recVecs = -Ng:Ng;
zeroVecs = recVecs * 0;

kap1n = kapn(q,om,g,Omega,e1,Ng);
kap2n = kapn(q,om,g,Omega,e2,Ng);
%produces the k vectors

out = zeros(2 * Ng + 1);
    for i = 1: 1
        %constructing vector of conductivity frequency coefficients
        sig = sigma2(om + recVecs * Omega, om, Omega, Ng)./ ...
              e0 .*(om + (recVecs).* Omega);
        %picks out an element in the conductivity vector to construct
        %diagonal
        sig_u = sig(Ng + 1 + i);%sigma2(om-i* Omega,om, Omega, Ng);
        sig_l = sig(Ng + 1 - i);%sigma2(om+i* Omega,om, Omega, Ng);
        
        %constructing diagonals
        sig_u = sig_u * (ones(1,2*Ng + 1 - i));
        sig_l = sig_l * (ones(1,2*Ng + 1 - i));
        upper_diag = 1i.*sig_u*alphaF; % upper diagonal;
        lower_diag = 1i.*sig_l*alphaF; %lower diagonal;
        out = out +  diag( upper_diag,+i) + diag( lower_diag,-i);
        %2N+1 dimension
    end
    %main diagonal
out = out + diag(e1./kap1n + e2./kap2n + 1i * sigma2(om + recVecs*Omega, om, Omega, Ng)./(e0*(om + recVecs*Omega)),0);

%  
% out = full(gallery('tridiag', ...
%             1i*sigmaDrude(om+(recVecs(2:end)-1)*Omega)*alphaF./ ...
%                 (e0*(om+(recVecs(2:end)-1)*Omega)) ,... % lower diagonal
%             (e1./kap1n+e2./kap2n) + 1i*sigmaDrude(om+recVecs*Omega)./ ...
%                 (e0*(om+recVecs*Omega)), ... % diagonal
%             1i*sigmaDrude(om+(recVecs(1:end-1)+1)*Omega)*alphaF./ ...
%                 (e0*(om+(recVecs(1:end-1)+1)*Omega)) )); % upper diagonal

end
%original
% out = full(gallery('tridiag', ...
%             1i*sigmaDrude(om+(recVecs(2:end)-1)*Omega)*alphaF./ ...
%                 (e0*(om+(recVecs(2:end)-1)*Omega)) ,... % lower diagonal
%             (e1./kap1n+e2./kap2n) + 1i*sigmaDrude(om+recVecs*Omega)./ ...
%                 (e0*(om+recVecs*Omega)), ... % diagonal
%             1i*sigmaDrude(om+(recVecs(1:end-1)+1)*Omega)*alphaF./ ...
%                 (e0*(om+(recVecs(1:end-1)+1)*Omega)) )); % upper diagonal



%modified code
% out = full(gallery('tridiag', ...
%             1i*sigmaDrude(om+(zeroVecs(2:end))*Omega)*alphaF./ ...
%                 (e0*(om+(zeroVecs(2:end))*Omega)) ,... % lower diagonal
%             (e1./kap1n+e2./kap2n) + 1i*sigmaDrude(om+zeroVecs*Omega)./ ...
%                 (e0*(om+zeroVecs*Omega)), ... % diagonal
%             1i*sigmaDrude(om+(zeroVecs(2:end))*Omega)*alphaF./ ...
%                 (e0*(om+(zeroVecs(2:end))*Omega)) )); % upper diagonal#

