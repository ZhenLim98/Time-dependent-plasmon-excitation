function [ out,kap1n,kap2n ] = trans_calc( q,om,g,Omega,e1,e2,Ng,alphaF,sou )

    [matrix_to_invert,kap1n,kap2n] = mMat(q,om,g,Omega,e1,e2,Ng,alphaF);
    out = matrix_to_invert\sou;

end

