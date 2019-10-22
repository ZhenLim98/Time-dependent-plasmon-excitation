function [ out,kap1n,kap2n ] = mMat( q,om,g,Omega,e1,e2,Ng,alphaF )
global e0

recVecs = -Ng:Ng;
zeroVecs = recVecs * 0;

kap1n = kapn(q,om,g,Omega,e1,Ng);
kap2n = kapn(q,om,g,Omega,e2,Ng);
%produces the k vectors

mod_conduc = sigmaDrude + alphaF / 2 * cos(Omega);

% out = zeros(2 * Ng + 1);
%     for i = 1:1
%         upper_diag = 1i*sigmaDrude(om+(recVecs((i + 1):end)-i)*Omega)*alphaF./ ...
%                      (e0*(om+(recVecs((i + 1):end)-i)*Omega));
%         lower_diag = 1i*sigmaDrude(om+(recVecs(1:(end-i))+i)*Omega)*alphaF./ ...
%                      (e0*(om+(recVecs(1:(end-i))+i)*Omega)); % upper diagonal
%         out = out +  diag( upper_diag,+i) + diag( lower_diag,-i);
%         %2N+1 dimension
%     end
% out = out + diag(e1./kap1n + e2./kap2n + 1i * sigmaDrude(om + recVecs*Omega)./(e0*(om + recVecs*Omega)),0);
 
out = full(gallery('tridiag', ...
            1i*sigmaDrude(om+(recVecs(2:end)-1)*Omega)*alphaF./ ...
                (e0*(om+(recVecs(2:end)-1)*Omega)) ,... % lower diagonal
            (e1./kap1n+e2./kap2n) + 1i*sigmaDrude(om+recVecs*Omega)./ ...
                (e0*(om+recVecs*Omega)), ... % diagonal
            1i*sigmaDrude(om+(recVecs(1:end-1)+1)*Omega)*alphaF./ ...
                (e0*(om+(recVecs(1:end-1)+1)*Omega)) )); % upper diagonal

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

