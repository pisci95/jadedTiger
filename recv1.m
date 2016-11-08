LL = 200; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 25; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma = 1; % Noise standard deviation. Default is 1

clc

%% message

p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ]; 
p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
p1 = p1.';

dataSize = 3;
timingSync  = [-1 1 -1 1 -1 1 1 -1 -1 1 1 1 -1 -1 -1 1 1 -1 -1 1 -1 1];
pilot = [-1, -1, -1, 1];
msgSize = (dataSize + length(timingSync) + length(pilot));
xtraSize = 2*msgSize;
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
[xcorred, lags] = xcorr(real(y1), real(toCorr));
%plot(xcorred);

[corrVal, offArg] = max(abs(xcorred));
offset = lags(offArg) + 125;
v1 = real(y1(offset:(offset+msgSize*fs)));


%% demodulate

w1 = flipud(p1); 
z1 = conv(w1, v1)*(1/fs);








%% sample it


%z1k = z1(ceil(Ns)+ceil(delay*fs):fs:end);
%z1k = z1(offset:fs:(msgSize*fs)+offset);
z1k = z1(125:fs:(msgSize*fs)+125);

%% guess


bits1_hat = sign(real(z1k));

res = bits1_hat(2:end);
correcting = res(msgSize-dataSize:end);
if (correcting(1) == 1)
    msg = correcting(2:end);
else
    msg = -1 .* correcting(2:end);
end

%% mapping

data = (msg + 1) .* .5



%% plot stuff

figure(1)
clf
subplot(4,2,1)
plot(real(y1),'b')
hold on
plot(imag(y1),'r')
legend('real','imag')
ylabel('yI(t)  and  yQ(t)')
xlabel('Time in samples')

subplot(4,2,2)
plot([0:length(y1)-1]/length(y1)-0.5, abs(fftshift(fft(y1))))
ylabel('abs(Y(f))')
xlabel('Frequency in 1/samples')

subplot(4,2,3)
plot(real(z1),'b')
hold on
plot(imag(z1),'r')
zoom xon
legend('real','imag')
ylabel('zI(t)  and  zQ(t)')
xlabel('Time in samples')


subplot(4,2,4)
plot([0:length(z1)-1]/length(z1)-0.5, abs(fftshift(fft(z1))))
ylabel('abs(Z(f))')
xlabel('Frequency in 1/samples')


subplot(4,2,5)
plot(real(z1k),'b')
hold on
plot(imag(z1k),'r')
zoom xon
legend('real','imag')
ylabel('z1kI(t)  and  z1kQ(t)')
xlabel('Time in samples')


subplot(4,2,6)
plot([0:length(z1k)-1]/length(z1k)-0.5, abs(fftshift(fft(z1k))))
ylabel('abs(Z1k(f))')
xlabel('Frequency in 1/samples')


subplot(4,2,7)
stem(res,'b')
%hold on
%stem(msg,'r')
zoom xon
%legend('real')
ylabel('bits(t)')
xlabel('Time in samples')

subplot(4,2,8)
plot([0:length(res)-1]/length(res)-0.5, abs(fftshift(fft(res))))
ylabel('abs(BITS(f))')
xlabel('Frequency in 1/samples')

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