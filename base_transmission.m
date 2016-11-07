
% Baseband modulation and demodulation matlab demo
% Two possible modulators demodulators, employing two different pulse shapes are shown

LL = 1000; % Total number of bits Default is 1000
T = 1; % Symbol period in microsec. Default is 1
N = 11; % length of filter in symbol periods. Default is 11
alpha = 0.2; % alpha of sqrt raised cosine filter (try 0.2, as well as 0.0 and 1.0)
nyquistflag = true; % true to use Nyquist pulse, false to use non-Nyquist pulse (stretch the Nyquist pulse)
fs = 100; % Over-sampling factor (Sampling frequency/symbol rate). Default is 100
Ns = floor(N*fs); % Number of filter samples
sigma_n = 1; % Noise standard deviation. Default is 1. Also try 5
delay = 0; % Offset from optimum sampling point (as fraction of symbol period). Try either 0 or 1/3

clc


% **********************************************************
% Non-Nyqusit pulse?
switch nyquistflag
    case true
        Tpulse_factor = 1.0; % Nyquist pulse
        disp('Using Nyquist pulses!')
    case false
        Tpulse_factor = 2.0; % False symbol period for pulse creates non-Nyquist pulse
        disp('WARNING! Using non-Nyquist pulses!')
        disp(' ')
end

% Initialize random number generator
rng(0);


% **********************************************************
% Modulation

% Create bits
bits = sign(randn(LL,1));

% Use rectangular pulse as one possible filter 
p1 = [zeros(ceil((Ns/2-fs/2)*Tpulse_factor),1) ; ones(ceil(fs*Tpulse_factor),1) ; zeros(Ns-fs-ceil((Ns/2-fs/2)*Tpulse_factor),1) ]; p1 = p1/norm(p1)/sqrt(1/(fs*Tpulse_factor)); % '1/fs' simply serves as 'delta' to approximate integral as sum

% Use sqrt-raised cosine filter form  ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as another possible filter
p2 = firrcos(Ns,1/2/T/Tpulse_factor,alpha,fs/T,'rolloff','sqrt'); p2 = p2/norm(p2)/sqrt(1/(fs*Tpulse_factor)); % '1/fs' simply serves as 'delta' to approximate integral as sum



% Create baseband signals
bit_up = upsample(bits,fs);
x1 = conv(bit_up,p1);
x2 = conv(bit_up,p2);
len = min([length(x1) length(x2)]);
x1 = x1(1:len); x2 = x2(1:len); 



% **********************************************************
% Channel

% Send them through separate AWGN channels
y1 = x1+sigma_n*randn(len,1); 
y2 = x2+sigma_n*randn(len,1); 



% **********************************************************
% Demodulation

% Matched filter
w1 = flipud(p1); 
w2 = flipud(p2);

% Filter with matched filter
z1 = conv(w1,y1)*(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum 
z2 = conv(w2,y2)*(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum 

% Sample filtered signal
z1k = z1(ceil(Ns)+ceil(delay*fs):fs:end); z1k = z1k(1:LL);
z2k = z2(ceil(Ns)+ceil(delay*fs):fs:end); z2k = z2k(1:LL);


% Detection
bits1_hat = sign(z1k);
bits2_hat = sign(z2k);


% Compute Bit error rate (BER)
BER1 = mean(bits1_hat ~= bits)
BER2 = mean(bits2_hat ~= bits)



% **********************************************************
% Alternative: Simple A/D converter-based detection (i.e., no matched filter)

% Sample y1, y2 directly (without matched filtering) and do detection
y1k = y1(ceil(Ns/2):fs:end); y1k = y1k(1:LL);
y2k = y2(ceil(Ns/2):fs:end); y2k = y2k(1:LL);
bits1_alt = sign(y1k);
bits2_alt = sign(y2k);

% Compute Bit error rate (BER) of simple A/D
BER1_alt = mean(bits1_alt ~= bits)
BER2_alt = mean(bits2_alt ~= bits)




% **********************************************************
% Waveforms and spectra

% Plot time domain signals
ax = [];
figure(1)
LargeFigure(gcf, 0.15); % Make figure large
clf
ax(1) = subplot(4,2,1);
plot([1:len]/fs*T,x1)
ylabel('$x_1(t)$')
axis tight
yl = ylim; ylim(yl*1.5);
ax(2) = subplot(4,2,2);
plot([1:len]/fs*T,x2)
ylabel('$x_2(t)$')
axis tight
ax(3) = subplot(4,2,3);
plot([1:len]/fs*T,y1)
ylabel('$y_1(t)$')
ax(4) = subplot(4,2,4);
plot([1:len]/fs*T,y2)
ylabel('$y_2(t)$')
ax(5) = subplot(4,2,5);
plot([1:length(z1)]/fs*T,z1)
ylabel('$z_1(t)$')
xlabel('time  in  microseconds')
ax(6) = subplot(4,2,6);
plot([1:length(z2)]/fs*T,z2)
ylabel('$z_2(t)$')
xlabel('time  in  microseconds')
ax(7) = subplot(4,2,7);
stem([1:LL],bits,'bx')
hold on
stem([1:LL],z1k,'ro')
ylabel('$x_k, z_{1,k}$')
xlabel('discrete time  $k$ (sampled at $t=kT$)')
axis tight
ax(8) = subplot(4,2,8);
stem([1:LL],bits,'bx');
hold on
stem([1:LL],z2k,'ro')
ylabel('$x_k, z_{2,k}$')
xlabel('discrete time  $k$ (sampled at $t=kT$)')
linkaxes(ax,'x')
axis tight
zoom xon

% Plot frequency domain signals
figure(2)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,2,1)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x1)))))
ylabel('$|X_1(f)|$')
axis([-4/T 4/T -40 40])
title('Frequency responses in dB')
subplot(3,2,2)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(x2)))))
ylabel('$|X_2(f)|$')
axis([-4/T 4/T -40 40])
title('Frequency responses in dB')
subplot(3,2,3)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(y1)))))
ylabel('$|Y_1(f)|$')
axis([-4/T 4/T -40 40])
subplot(3,2,4)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(y2)))))
ylabel('$|Y_2(f)|$')
axis([-4/T 4/T -40 40])
subplot(3,2,5)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(1/sqrt(len)*fft(z1(1:len))))))
ylabel('$|Z_1(f)|$')
axis([-4/T 4/T -40 40])
xlabel('frequency  in  MHz')
subplot(3,2,6)
plot([-len/2+1:len/2]/len*fs/T,20*log10(abs(fftshift(fft(1/sqrt(len)*z2(1:len))))))
ylabel('$|Z_2(f)|$')
axis([-4/T 4/T -40 40])
xlabel('frequency  in  MHz')
figure(1)
