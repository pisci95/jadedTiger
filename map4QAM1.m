function [ mapped ] = map4QAM1( input, d )
%4QAM1 maps the 4QAM values
%   Detailed explanation goes here

input = [input, zeros(1, mod(length(input),2))];
lenOut = (length(input)/2);
mapped = zeros(1, lenOut);
for i = 1:lenOut
    if (input(i:i+1) == [0 0])
        mapped(i) = -1 - 1i;
    elseif (input(i:i+1) == [0 1])
        mapped(i) = -1 + 1i;
    elseif (input(i:i+1) == [1 0])
        mapped(i) = 1 - 1i;
    else 
        mapped(i) = 1 + 1i;
    end
end
mapped = mapped .* (d/2);


end

