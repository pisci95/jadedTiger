function [msgsym] = OFDMDemod(numSubcarriers, D, symbols)
%% OFDM Demodulator: Recovers signals from OFDM
% [numSubcarriers] subcarriers and [D] cyclic prefix length

% FFT
subc = numSubcarriers; % #subcarriers
fftmsg = zeros(ceil(length(symbols)/(subc+D)), subc);

% take fft of each subc+D-sized chunk
for ii = 1:ceil(length(symbols)/(subc+D))
    % Discard first D samples
    fftmsg(ii, 1:subc) = fft(symbols((subc+D)*(ii-1)+1+D:(subc+D)*(ii-1)+subc+D), subc); 
end
% last = symbols((subc+D)*(ii)+D+1:end);
% % zero pad last entry in case of unidealities
% fftmsg(ii+1, 1:subc) = fft(last, subc);

width = length(fftmsg(:, 1));
msgsym = zeros(1, width * subc);
% concatenate them all into a message
for ii = 1:width
    msgsym((ii-1)*subc+1:ii*subc) = fftmsg(ii,:);
end

end