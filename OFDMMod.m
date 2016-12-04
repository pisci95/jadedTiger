function [toTransmit] = OFDMMod(numSubcarriers, D, msgSym)
% 
% %% OFDM Modulator: modulates msgSym into OFDM subcarriers
% % [numSubcarriers] subcarriers and [D] cyclic prefix length
% 
% % IFFT
% subc = numSubcarriers; % #subcarriers
% ifftmsg = zeros(ceil(length(msgSym)/subc), subc + D);
% 
% % break symbols into chunks of size numSubcarriers
% for ii = 1:floor(length(msgSym)/subc)
%     % take fft of each chunk and add cyclic prefix
%     ifftmsg(ii, D+1:D+subc) = ifft(msgSym(subc*(ii-1)+1:subc*(ii-1)+subc), subc); 
%     ifftmsg(ii, 1:D) = ifftmsg(ii, end-D+1:end);
%     fftmsg(ii, 1:subc) = fft(ifftmsg(ii, D+1:subc+D), subc); 
% end
% last = msgSym(subc*(ii)+1:end);
% % zero pad last entry, do the same
% ifftmsg(ii+1, D+1:D+subc) = ifft([last zeros(1, subc-length(last))], subc);
% ifftmsg(ii+1, 1:D) = ifftmsg(ii+1, end-D+1:end);
% fftmsg(ii+1, 1:subc) = fft(ifftmsg(ii+1,D+1:end), subc);
% 
% width = length(fftmsg(:, 1));
% msgsym = zeros(1, width * subc);
% % concatenate them all into a message
% for ii = 1:width
%     msgsym((ii-1)*subc+1:ii*subc) = fftmsg(ii,:);
% end
% 
% end

%% OFDM Modulator: modulates msgSym into OFDM subcarriers
% [numSubcarriers] subcarriers and [D] cyclic prefix length

% IFFT
subc = numSubcarriers; % #subcarriers
ifftmsg = zeros(ceil(length(msgSym)/subc), subc + D);

% break symbols into chunks of size numSubcarriers
for ii = 1:floor(length(msgSym)/subc)
    % take fft of each chunk and add cyclic prefix
    ifftmsg(ii, D+1:D+subc) = ifft(msgSym(subc*(ii-1)+1:subc*(ii-1)+subc), subc); 
    ifftmsg(ii, 1:D) = ifftmsg(ii, end-D+1:end);
end
last = msgSym(subc*(ii)+1:end);
% zero pad last entry, do the same
ifftmsg(ii+1, D+1:D+subc) = ifft([last zeros(1, subc-length(last))], subc);
ifftmsg(ii+1, 1:D) = ifftmsg(ii+1, end-D+1:end);

width = length(ifftmsg(:, 1));
toTransmit = zeros(1, width * (subc+D));
%concatenate them all into a time signal
for ii = 1:width
    toTransmit((ii-1)*(subc+D)+1:ii*(subc+D)) = ifftmsg(ii,:);
end