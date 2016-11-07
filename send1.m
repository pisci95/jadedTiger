
LL = 200; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma = 1; % Noise standard deviation. Default is 1

clc

%% message
%bits = sign(randn(1, LL));

dataSize = 3;
bits = [1 1 -1];

timingSync  = [1 -1 1 -1 1 1 -1 -1 1 1 1 -1 -1 -1 1 1 -1 -1 1 -1 1];
pilot = [-1, -1, 1, -1];
buffer = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];

msg = [buffer, timingSync, pilot, bits];



%% impulse
p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ]; 
p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
p1 = p1.';

%% Create baseband signals
bit_up = upsample(msg,fs);
x1 = conv(bit_up,p1);
x2 = conv(bit_up,p1);
x3 = x1 + 1i .* x2;
transmitsignal = x1.';
%len = length(x1);
%plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x1)))))

save('transmitsignal.mat', 'transmitsignal');