function [ decoded ] = guess4QAM1( pkt )
%guess4QAM1 Summary of this function goes here
%   Detailed explanation goes here

lenOut = (length(pkt));
decoded = zeros(1, 2*lenOut);
for i = 1:lenOut
    valR = real(pkt(i));
    valI = imag(pkt(i));
    
    if ((valR < 0) && (valI < 0))
        decoded(2*i - 1) = 0;
        decoded(2*i) = 0;
    elseif ((valR < 0) && (valI > 0))
        decoded(2*i - 1) = 0;
        decoded(2*i) = 1;
    elseif ((valR > 0) && (valI < 0))
        decoded(2*i - 1) = 1;
        decoded(2*i) = 0;
    else 
        decoded(2*i - 1) = 1;
        decoded(2*i) = 1;
    end
end

end

