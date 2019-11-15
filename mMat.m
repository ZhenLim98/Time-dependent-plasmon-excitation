function [ out,kap1n,kap2n ] = mMat( q,om,g,Omega,e1,e2,Ng )
global e0

cl = 3e8;
alphaF = 1;
%increasing the strength of modulation
recVecs = -Ng:Ng;

kap1n = kapn(q,om,g,Omega,e1,Ng);
kap2n = kapn(q,om,g,Omega,e2,Ng);
%produces the k vectors
out = zeros(2 * Ng + 1);
    for i = 1:Ng
        %sigma2 produces a row vector of conductivity coefficients that are
        %also a function of frequency
        upper_vec = lambda_i(recVecs, recVecs((i + 1):end), -i, Ng, Omega)*alphaF ;
        %the same all the way across the diagonal
        upper_diag = -cl/(om).*upper_vec./ ...
                     (e0*(om+(recVecs((i + 1):end)-i)*Omega));
        lower_vec = lambda_i(recVecs, recVecs(1:(end-i)),+i, Ng, Omega) *alphaF ;
        %the same all the way across the diagonal
        lower_diag = -cl/(om).*lower_vec./ ...
                     (e0*(om+(recVecs(1:(end-i))+i)*Omega)); % upper diagonal   
        out = out +  diag( upper_diag,+i) + diag( lower_diag,-i);
        %2N+1 dimension
    end
    
out = out + diag(e1./kap1n + e2./kap2n - cl/om *alphaF* lambda_i(recVecs, recVecs, 0, Ng, Omega)./(e0*(om + recVecs*Omega)),0);

end