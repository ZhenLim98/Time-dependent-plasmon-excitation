function [ out,kap1n,kap2n ] = mMat2( q,om,g,Omega,e1,e2,Ng,alphaF )
global e0 cl

recVecs = -Ng:Ng;
zeroVecs = recVecs * 0;

kap1n = kapn(q,om,g,Omega,e1,Ng);
kap2n = kapn(q,om,g,Omega,e2,Ng);
%produces the k vectors
out = zeros(2 * Ng + 1);
    for i = 1:Ng
        %sigma2 produces a row vector of conductivity coefficients that are
        %also a function of frequency
        upper_vec = lambda( (om+(recVecs((i + 1):end)-i)*Omega),om, Omega,(recVecs((i + 1):end)), Ng);
        upper_diag =  - cl/(om*e0) .* upper_vec * alphaF./ ...
                     (e0*(om+(recVecs((i + 1):end)-i)*Omega));
        lower_vec = lambda(om+(recVecs(1:(end-i))+i)*Omega, om,Omega, (recVecs(1:(end-i))), Ng);
        lower_diag = - cl/(om*e0) .* lower_vec * alphaF./ ...
                     (e0*(om+(recVecs(1:(end-i))+i)*Omega)); % upper diagonal
        out = out +  diag( upper_diag,+i) + diag(lower_diag,-i);
        %2N+1 dimension
    end
out = out + diag(e1./kap1n + e2./kap2n - cl * lambda(om + recVecs*Omega, om, Omega, recVecs, Ng )./(om* e0 *(om + zeroVecs*Omega)),0);
%Check upper lines
end
