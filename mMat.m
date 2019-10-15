function [ out,kap1n,kap2n ] = mMat( q,om,g,Omega,e1,e2,Ng,alphaF )
global e0

recVecs = -Ng:Ng;

kap1n = kapn(q,om,g,Omega,e1,Ng);
kap2n = kapn(q,om,g,Omega,e2,Ng);

out = full(gallery('tridiag', ...
            1i*sigmaDrude(om+(recVecs(2:end)-1)*Omega)*alphaF./ ...
                (e0*(om+(recVecs(2:end)-1)*Omega)) ,... % lower diagonal
            (e1./kap1n+e2./kap2n) + 1i*sigmaDrude(om+recVecs*Omega)./ ...
                (e0*(om+recVecs*Omega)), ... % diagonal
            1i*sigmaDrude(om+(recVecs(1:end-1)+1)*Omega)*alphaF./ ...
                (e0*(om+(recVecs(1:end-1)+1)*Omega)) )); % upper diagonal

end

