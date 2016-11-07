% modulation matlab demo
% Two possible modulators, employing two different pulse shapes are shown

LL = 1000; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0; % alpha of sqrt raised cosine filter
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma_n = 1; % Noise standard deviation. Default is 1

clc


% Create bits
bits = sign(randn(LL,1));

% Use rectangular pulse as one possible filter 
p1 = [zeros(ceil(Ns/2-fs/2),1) ; ones(fs,1) ; zeros(Ns-fs-ceil(Ns/2-fs/2),1) ];
p1 = p1/norm(p1)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum

% Use sqrt-raised cosine filter form  ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as another possible filter
p2 = firrcos(Ns,1/2/T,alpha,fs/T,'rolloff','sqrt');
p2 = p2/norm(p2)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum



% Create baseband signals
bit_up = upsample(bits,fs);
x1 = conv(bit_up,p1);
x2 = conv(bit_up,p2);
len = min([length(x1) length(x2)]);
x1 = x1(1:len); x2 = x2(1:len); 


% Send them through separate AWGN channels
y1 = x1+sigma_n*randn(len,1); 
y2 = x2+sigma_n*randn(len,1); 


% Plot time domain signals
close all
ax = [];
figure(1)
%LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,2,1)
plot([-Ns/2+1:Ns/2]/fs*T,p1)
ylabel('$p_1(t)$')
axis auto
yl = ylim; ylim(yl*1.5);
subplot(3,2,2)
plot([-Ns/2:Ns/2]/fs*T,p2)
ylabel('$p_2(t)$')
axis tight
ax(1) = subplot(3,2,3);
plot([1:len]/fs*T,x1)
ylabel('$x_1(t)$')
axis tight
yl = ylim; ylim(yl*1.5);
ax(2) = subplot(3,2,4);
plot([1:len]/fs*T,x2)
ylabel('$x_2(t)$')
axis tight
ax(3) = subplot(3,2,5);
plot([1:len]/fs*T,y1)
ylabel('$y_1(t)$')
xlabel('time  in  microseconds')
axis tight
ax(4) = subplot(3,2,6);
plot([1:len]/fs*T,y2)
ylabel('$y_2(t)$')
xlabel('time  in  microseconds')
axis tight
linkaxes(ax,'x')
zoom xon


% Plot frequency domain signals
figure(2)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,2,1)
plot([-Ns/2+1:Ns/2]/Ns*fs/T,20*log10(abs(fftshift(1/sqrt(Ns)*fft(p1)))))
ylabel('$|P_1(f)|$')
axis([-4/T 4/T -40 40])
title('Frequency responses in dB')
subplot(3,2,2)
plot([-Ns/2:Ns/2]/Ns*fs/T,20*log10(abs(fftshift(1/sqrt(Ns)*fft(p2)))))
ylabel('$|P_2(f)|$')
axis([-4/T 4/T -40 40])
title('Frequency responses in dB')
subplot(3,2,3)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x1)))))
ylabel('$|X_1(f)|$')
axis([-4/T 4/T -40 40])
subplot(3,2,4)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x2)))))
ylabel('$|X_2(f)|$')
axis([-4/T 4/T -40 40])
subplot(3,2,5)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(y1)))))
ylabel('$|Y_1(f)|$')
axis([-4/T 4/T -40 40])
xlabel('frequency  in  MHz')
subplot(3,2,6)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(y2)))))
ylabel('$|Y_2(f)|$')
axis([-4/T 4/T -40 40])
xlabel('frequency  in  MHz')
figure(1)
