
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
bits = [0 1 0];

timingSync  = [0 1 0 1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 0 1 0 1];
pilot = [0, 0, 0, 1];
buffer = [ones(1, 300), zeros(1, 300)];

msg = [timingSync, pilot, bits, pilot];


%% symbol mapping
msgSym = msg .* 2 -1;
msgSym = [buffer -1 -1 -1 msgSym] .* .75;



%% impulse
p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ]; 
p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
p1 = p1.';

%% Create baseband signals
bit_upR = upsample(real(msgSym),fs);
bit_upI = upsample(imag(msgSym),fs);
xR = conv(bit_upR,p1);
xI = conv(bit_upI,p1);
x1 = xR + 1i .* xI;
transmitsignal = x1.';
%len = length(x1);
%plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x1)))))


%% plot stuff

figure(1)
clf
subplot(4,2,1)
plot(p1,'b')
%hold on
%plot(imag(y1),'r')
%legend('real','imag')
ylabel('p(t)')
xlabel('Time in samples')

subplot(4,2,2)
plot([0:length(p1)-1]/length(p1)-0.5, abs(fftshift(fft(p1))))
ylabel('abs(P(f))')
xlabel('Frequency in 1/samples')

subplot(4,2,3)
stem(msg,'b')
zoom xon
ylabel('msgk')
xlabel('signal in bits')


subplot(4,2,4)
plot([0:length(msg)-1]/length(msg)-0.5, abs(fftshift(fft(msg))))
ylabel('abs(MSG(f))')
xlabel('Frequency in 1/samples')


subplot(4,2,5)
plot(real(msgSym),'b')
hold on
plot(imag(msgSym),'r')
zoom xon
legend('real','imag')
ylabel('mapI(t) & mapQ(t)')
xlabel('Time in samples')


subplot(4,2,6)
plot([0:length(msgSym)-1]/length(msgSym)-0.5, abs(fftshift(fft(msgSym))))
ylabel('abs(MAP(f))')
xlabel('Frequency in 1/samples')


subplot(4,2,7)
plot(real(x1),'b')
hold on
plot(imag(x1),'r')
zoom xon
legend('real' ,'imag')
ylabel('xI(t) & xQ(t)')
xlabel('Time in samples')

subplot(4,2,8)
plot([0:length(x1)-1]/length(x1)-0.5, abs(fftshift(fft(x1))))
ylabel('abs(X(f))')
xlabel('Frequency in 1/samples')


save('transmitsignal.mat', 'transmitsignal');