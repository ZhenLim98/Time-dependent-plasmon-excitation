function [ ] = check_subwavelength(g,k0)

if k0 > g
    warning('Grating is not subwavelength!')
end

