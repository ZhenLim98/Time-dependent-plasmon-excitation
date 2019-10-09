function [k1 , k2] = getk_2(n, G, q, e_1, e_2, omega, c)

G_vector = ( -n : n ) * G;
k1 =  -sqrt( (q +  G_vector) .^ 2 - e_1 * omega * omega / (c * c) );
k2 = -1i * abs(sqrt( (q +  G_vector) .^ 2 - e_2 * omega * omega / (c * c) ));

end