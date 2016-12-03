function [encodedBits] = encode(message)
state = '00';
for ii = 1:length(message)
    bit = message(ii);
    switch(state)
        case '00'
            if bit == 1
                encodedBits((2*ii - 1):2*ii) = [1 1];
                state = '01';
            else % bit is 0
                encodedBits((2*ii - 1):2*ii) = [0 0];
                state = '00';
            end
        case '01'
            if bit == 1
                encodedBits((2*ii - 1):2*ii) = [0 1];
                state = '11';
            else % bit is 0
                encodedBits((2*ii - 1):2*ii) = [1 0];
                state = '10';
            end
        case '10'
            if bit == 1
                encodedBits((2*ii - 1):2*ii) = [0 0];
                state = '01';
            else % bit is 0
                encodedBits((2*ii - 1):2*ii) = [1 1];
                state = '00';
            end
        case '11'
            if bit == 1
                encodedBits((2*ii - 1):2*ii) = [1 0];
                state = '11';
            else % bit is 0
                encodedBits((2*ii - 1):2*ii) = [0 1];
                state = '10';
            end
        otherwise
            %cry %should not reach
    end
end