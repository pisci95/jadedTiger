% Viterbi Decoder
function [message] = decode(receivedBits)

message = zeros(1, length(receivedBits)/2);

%% Step 1: Allocate Table of 4x[bits]/rate weights
weights = zeros(8, length(receivedBits)/2);
% rows 1-8: 00-00, 00-01, 01-10, 01-11, 10-00, 10-01, 11-10, 11-11
transition = zeros(8, 2);
transition(1, :) = [0 0];
transition(2, :) = [1 1];
transition(3, :) = [1 0];
transition(4, :) = [0 1];
transition(5, :) = [1 1];
transition(6, :) = [0 0];
transition(7, :) = [0 1];
transition(8, :) = [1 0];

% First transitions: must be state 00
weights(1, 1) = abs(receivedBits(1)-0) + abs(receivedBits(2)-0);
weights(2, 1) = abs(receivedBits(1)-1) + abs(receivedBits(2)-1);
weights(3:8, 1) = inf;

% For each two bits
for ii = 2:length(receivedBits)/2
    % Calculate Hamming distance for each transition and store in array
    bits = receivedBits(2*ii-1: 2*ii);
    for jj = 1:8
        dist = abs(bits(1)-transition(jj, 1)) + abs(bits(2)-transition(jj, 2));
        weights(jj, ii) = dist;
    end
end

%% Step 2: Calculate path lengths
paths = zeros(4, length(receivedBits)/2);

% initial path lengths
paths(1, 1) = weights(1, 1);
paths(2, 1) = weights(2, 1);
paths(3:4, 1) = inf;

% for each two bits
for ii = 2:length(receivedBits)/2
    % calculate path lengths for this instant
    paths(1, ii) = min((paths(1, ii-1) + weights(1, ii)),(paths(3, ii-1) + weights(5, ii)));
    paths(2, ii) = min((paths(1, ii-1) + weights(2, ii)),(paths(3, ii-1) + weights(6, ii)));
    paths(3, ii) = min((paths(2, ii-1) + weights(3, ii)),(paths(4, ii-1) + weights(7, ii)));
    paths(4, ii) = min((paths(2, ii-1) + weights(4, ii)),(paths(4, ii-1) + weights(8, ii)));
end

% find most likely final state
[minimum, index] = min([paths(1, length(receivedBits)/2)
paths(2, length(receivedBits)/2)
paths(3, length(receivedBits)/2)
paths(4, length(receivedBits)/2)]);

state = [0 0];

if index == 1
    state = [0 0];
elseif index == 2
    state = [0 1];
elseif index == 3
    state = [1 0];
elseif index == 4
    state = [1 1];
end

% Trace back to get bits
for ii = length(receivedBits)/2:-1:1
    if isequal(state, [0 0])
        message(ii) = 0;
        if ii > 1
            if (paths(1, ii-1) + weights(1, ii) == paths(1, ii))
                state = [0 0];
            elseif (paths(3, ii-1) + weights(5, ii) == paths(1, ii))
                state = [1 0];
            end
        end
    elseif isequal(state, [0 1])
        message(ii) = 1;
        if ii>1
            if (paths(1, ii-1) + weights(2, ii) == paths(2, ii))
                state = [0 0];
            elseif (paths(3, ii-1) + weights(6, ii) == paths(2, ii))
                state = [1 0];
            end
        end
    elseif isequal(state, [1 0])
        message(ii) = 0;
        if ii>1
            if (paths(2, ii-1) + weights(3, ii) == paths(3, ii))
                state = [0 1];
            elseif (paths(4, ii-1) + weights(7, ii) == paths(3, ii))
                state = [1 1];
            end
        end
    elseif isequal(state, [1 1])
        message(ii) = 1;
        if ii>1
            if (paths(2, ii-1) + weights(4, ii) == paths(4, ii))
                state = [0 1];
            elseif (paths(4, ii-1) + weights(8, ii) == paths(4, ii))
                state = [1 1];
            end
        end
    end
end