
LL = 10; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma_n = 1; % Noise standard deviation. Default is 1

clc

%% message
bits = sign(randn(LL,1));

%% impulse
p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ]; p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum


%% Create baseband signals
bit_up = upsample(bits,fs);
x1 = conv(bit_up,p1);



%% Send them through separate AWGN channels
a = randn(size(x1));
b = randn(size(x1));
noiseZ= sigma/sqrt(2)*(a + j*b);

y1 = noiseZ + x1;






%% plot stuff