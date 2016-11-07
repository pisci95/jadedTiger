
LL = 200; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma = 1; % Noise standard deviation. Default is 1

clc

%% message
bits = sign(randn(1, LL));


timingSync  = [1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1 ];
pilot = [-1, -1, -1,-1, -1, -1,-1, -1, -1, 1];

msg = [timingSync, pilot, bits];


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

%% Send them through separate AWGN channels


%load receivedsignal.mat;

% Received signal
if exist('receivedsignal.mat','file')
    load('receivedsignal.mat');
    y1 = receivedsignal.';
else
    disp(['Error! Did not find receivedsignal.mat file.'])
    return
end

if ~exist('receivedsignal','var')
    disp('Error! Loaded file does not contain the receivedsignal variable.')
    return;
end
% a = randn(size(x1));
% b = randn(size(x1));
% noiseZ= sigma/sqrt(2)*(a + j*b);
% 
% y1 = noiseZ + x1;



%% timing recovery

time_up = upsample(timingSync,fs);
toCorr = conv(time_up, p1);
[xcorred, lags] = xcorr(y1, toCorr);
%plot(xcorred);

[corrVal, offArg] = max(abs(xcorred));
delay = lags(offArg);



%% demodulate

w1 = flipud(p1); 
z1 = conv(w1, y1)*(1/fs);


%% sample it


z1k = z1(ceil(Ns)+ceil(delay*fs):fs:end);
z1k = z1k(1:225);


%% guess


bits1_hat = sign(real(z1k));


%% BER


BER1 = mean(bits1_hat ~= msg)

%% plot stuff

figure(1)
clf
subplot(4,1,1)
plot(real(x1),'b')
hold on
plot(imag(x1),'r')
legend('real','imag')
ylabel('xI(t)  and  xQ(t)')
xlabel('Time in samples')
subplot(4,1,2)
plot(real(y1),'b')
hold on
plot(imag(y1),'r')
zoom xon
legend('real','imag')
ylabel('yI(t)  and  yQ(t)')
xlabel('Time in samples')

subplot(4,1,3)
plot(real(z1),'b')
hold on
plot(imag(z1),'r')
zoom xon
legend('real','imag')
ylabel('zI(t)  and  zQ(t)')
xlabel('Time in samples')

subplot(4,1,4)
plot(lags, real(xcorred),'b')
hold on
plot(lags, imag(xcorred),'r')
zoom xon
legend('real','imag')
ylabel('xcorredI(t)  and  xcorredQ(t)')
xlabel('Time in samples')

% figure(2)
% clf
% subplot(2,1,1)
% plot([0:length(x1)-1]/length(x1)-0.5, abs(fftshift(fft(x1))))
% ylabel('abs(X(f))')
% xlabel('Frequency in 1/samples')
% subplot(2,1,2)
% plot([0:length(y1)-1]/length(y1)-0.5, abs(fftshift(fft(y1))))
% ylabel('abs(Y(f))')
% xlabel('Frequency in 1/samples')