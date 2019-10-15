function [ out ] = sigmaDrude( om )
global ee eF hbar gamma

out = ee^2/pi/hbar^2*eF./(gamma-1i*om);

end

