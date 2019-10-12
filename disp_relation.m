function [q_light, qdata, wdata ] = disp_relation(energy, spacing, w_bot,w_top, gamma)
%energy is in 
h_bar = 6.63e-34/(2 * pi);
c = 3e8;
w_p = energy/h_bar;
e_1 = 1;

nSteps = round((w_top - w_bot)/spacing);
wdata = linspace(0, 0, nSteps);
qdata = linspace(0, 0, nSteps);

    %constructing graph
    for frequency = w_bot : spacing : w_top

    w = 2 * pi * frequency;

    index = (frequency -  w_bot )/spacing + 1;

    e_2 =  1 - w_p * w_p/( w * w +  1i * w * gamma );

    qdata(index) = w/c * sqrt(e_1 * e_2/(e_1 + e_2));
    wdata(index) = w;

    end

%constructing lightline
    q_light = wdata / c;

end