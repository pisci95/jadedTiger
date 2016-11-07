
LL = 200; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma = 1; % Noise standard deviation. Default is 1

clc

%% message
bits = ones(1, LL);%sign(randn(1, LL));


timingSync  = [1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1 ];
pilot = [-1, -1, -1,-1, -1, -1,-1, -1, -1, 1];

msg = bits;%[timingSync, pilot, bits];


%% impulse
p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ]; 
p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
p1 = p1.';
p2 = p1;

%% Create baseband signals
bit_up = upsample(msg,fs);
x1 = conv(bit_up,p1);
transmitsignal = x1.';

save('transmitsignal.mat', 'transmitsignal');