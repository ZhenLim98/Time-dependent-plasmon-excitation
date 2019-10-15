function [ out ] = kapn( q,om,g,Omega,eps,Ng )
global cl

recVecs = -Ng:Ng;

out =  sqrt((q+recVecs*g).^2-eps*(om+recVecs*Omega).^2/cl^2) ...
            .* ( ((q+recVecs*g).^2) > (eps*(om+recVecs*Omega).^2/cl^2) ) ...
            - 1i*sqrt(eps*(om+recVecs*Omega).^2/cl^2-(q+recVecs*g).^2) ...
            .* ( ((q+recVecs*g).^2) < (eps*(om+recVecs*Omega).^2/cl^2) );

end

