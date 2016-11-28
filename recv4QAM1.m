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


pilotLen = 8;
msgLen = 32;
pktBodySz = 20;
d= 1;
pilot = (d/2)*(ones(1, pilotLen) + 1i*zeros(1, pilotLen)); 
timingSync  = (d/2)*(textEncoder('timeSync')*2 -1)

dataSize = 3;
msgSize = (dataSize + length(timingSync) + length(pilot));
xtraSize = 2*msgSize;
%load receivedsignal.mat;

% Received signal
if exist('receivedsignal.mat','file')
    load('receivedsignal.mat');
    y1 = receivedsignal.';
elseif exist('transmitsignal.mat','file')
    load('transmitsignal.mat');
    y1 = transmitsignal.';
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
%v1 = real(y1(offset:(offset+msgSize*fs)));
v1 = y1(offset:end); %new one
v2 = v1((25+length(timingSync))*fs:end);




%% plot received bits
% figure(1)
% hold on;
% title('Signal Space Diagram');
% xlabel('I');
% ylabel('Q');
% scatter(real(y1(offset:end)), imag(y1(offset:end)));


%% demodulate

w1 = transpose(conj(p1));
z1 = conv(w1, v1)*(1/fs);



%% sample it


%z1k = z1(ceil(Ns)+ceil(delay*fs):fs:end);
%z1k = z1(offset:fs:(msgSize*fs)+offset);
%z1k = z1(125:fs:(msgSize*fs)+125);
z1k = z1(125:fs:end);
v1k = v1(1:fs:end);

%% guess



%% plot sampled bits
%scatter(real(z1k), imag(z1k), 'rx');
%legend('y1', 'z1k');


ppSectZ = z1k(length(timingSync)+1:end);
ppSectV = v1k(length(timingSync)+1:end);



pktNum = ceil(5*msgLen/pktBodySz);
body = zeros(1, pktNum*(pktBodySz));
mpdSize = pktBodySz/2;
for i = 1:pktNum
    i
    bdyOff = (mpdSize)*(i-1);
    msgOff = (mpdSize+pilotLen)*(i-1);
    toEqZ = ppSectZ(msgOff+1:mpdSize+pilotLen+msgOff);
    toEqV = ppSectV(msgOff+1:mpdSize+pilotLen+msgOff);
    
    clf
    figure(5)
    subplot(4,1,1)
    plot(real(toEqZ),'b')
    hold on
    plot(imag(toEqZ),'r')
    %legend('real','imag')
    ylabel('toEqZ  and  bodyQ')
    xlabel('Time in samples')
    
    subplot(4,1,2)
    plot(real(toEqV),'b')
    hold on
    plot(imag(toEqV),'r')
    %legend('real','imag')
    ylabel('toEqZ  and  bodyQ')
    xlabel('Time in samples')
    
    eqed = equalize1(toEqZ, toEqV, pilotLen);
    body(bdyOff + 1: bdyOff + mpdSize) = eqed;
    subplot(4,1,3)
    plot(real(eqed),'b')
    hold on
    plot(imag(eqed),'r')
    %legend('real','imag')
    ylabel('toEqZ  and  bodyQ')
    xlabel('Time in samples')
    subplot(4,1,4)
    plot(real(ppSectZ),'b')
    hold on
    plot(imag(ppSectZ),'r')
    %legend('real','imag')
    ylabel('toEqZ  and  bodyQ')
    xlabel('Time in samples')
    
end

%plot z1k
figure(1)
plot(real(body),'b')
hold on
plot(imag(body),'r')
%legend('real','imag')
ylabel('body  and  bodyQ')
xlabel('Time in samples')

%bits1_hat = sign(real(z1k));
bits1_hat = guess4QAM1(z1k);


res = bits1_hat(2:end);
correcting = res(msgSize-dataSize:end);
if (correcting(1) == 1)
    msg = correcting(2:end);
else
    msg = -1 .* correcting(2:end);
end

%% mapping

data = (msg + 1) .* .5;



%% plot stuff

figure(3)
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
plot(real(ppSect),'b')
hold on
plot(imag(ppSect),'r')
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