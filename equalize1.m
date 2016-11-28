function [ v ] = equalize1( z, x, pilotLen )
%equalize1 takes in z, x, and pilot length and does 1-tap equalizing
%   Detailed explanation goes here


zp = z(1:pilotLen);
zr = z((pilotLen + 1):end);
xp = x(1:pilotLen);
temp= zp*(transpose(conj(xp)));
h =  temp /(norm(xp)^2);

v = zr ./ h;



end

