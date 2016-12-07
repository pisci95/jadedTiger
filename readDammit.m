% Process received signal to decode message
% Mid-project solution
% Only decodes one bit
%
% Hints for the final part of the project are provided as ** hints **
% Be sure to look closely at these hints, else you will have a tough time in the final part!
% It is good coding practice to define functions (e.g., bits2symbolMap(), pulseShaping(), etc.) to organize the code logically


clc
%close all
set(0,'DefaultTextInterpreter','latex')
disp(' ')

% ************************************************
% Fixed (given) parameters. THESE MUST NOT BE CHANGED
RXusersamplingrate = 25e6; % user TX baseband sampling rate in Hz. 


%% ************************************************
% RUNTIME or TEST?
flag = false; % true indicates normal operation. false indicates the ideal case of receivedsignal being substituted for transmit signal directly (to help debugging your code)

if flag == true
    % ************************************************
    % Load latest (based on time stamp) received file
    disp('Loading the received signal')
    fname = dir('receivedsignal*');
    fnamereceived = fname(end).name;
    load(fnamereceived) 
else
    % ************************************************
    % Testing code only...
    disp('CAUTION! Using a SIMULATED channel.')    

    transmitsignal = awgn(transmitsignal, 2);
    downsamplerate = ceil(TXusersamplingrate/RXusersamplingrate);
    receivedsignal = [zeros(1000,1); transmitsignal(1:downsamplerate:end); zeros(1000,1)]; % No noise. Downsampled TX signal used so as to received at the RX sampling rate

end


rxsig = reshape(receivedsignal,[],1); % Received complex signal yI(t)+jyQ(t)
disp(' ')



%% ************************************************
% Other parameters



%% ************************************************
% Upsample received signal to convert to higher sampling rate. 
upsamplingfactor = ceil(TXusersamplingrate/RXusersamplingrate); % Factor by which to upsample received signal. Chosen as the TX sampling rate for convenience. Finer synchronization will be achieved if a higher upsampled rate is used. ceil() does not matter due to choice of RX and TX rates
rxsigupsampled = resample(rxsig, RXusersamplingrate*upsamplingfactor, RXusersamplingrate); % Upsampled received signal




%% ************************************************
% Timing recovery as well as frame synchronization

CTtau = conv(rxsigupsampled,flipud(conj(syncSignal)));
[mxsync, indx] = max(abs(CTtau));
 % This marks the beginning of the sync pilot in the received signal. 


e = 5:1:1.5*(lenSyncSpl + lenBufSpl);% length(rxsigupsampled)-lenSyncSpl-1;
res =zeros(length(e), 1);
c = 1;
for i = e
    possync = i;
    tp = abs(corrcoef(rxsigupsampled(possync:possync+lenSyncSpl-1),syncSignal));
    synccorrcoef = tp(2,1);
    res(c) = synccorrcoef;
    c = c+1;
end
[m2, i2] = max(abs(res));
possync = e(i2)


tp = abs(corrcoef(rxsigupsampled(possync:possync+lenSyncSpl-1),syncSignal));
synccorrcoef = tp(2,1);
disp(['Correlation coefficient of synchronization is  ' num2str(synccorrcoef)])
if synccorrcoef > 0.8
    disp('So, sync seems to have worked well!')
else
    disp('Check sync again. May not have worked well. I will proceed anyway!')
end




%% ************************************************
% Time-offset received signal, based on extracted timing.
% Also does frame synchronization correctly, since the pilot bits' autocorrelation
% function is sharp


%removes buffer
yt = rxsigupsampled(possync:possync+lenTXSpl-1);



%% ************************************************
% Extract signal outside our own signal interval

nonsignalpre  = rxsigupsampled(1:possync-1); % Waveform prior to our signal interval
nonsignalpost  = rxsigupsampled(possync+lenTXSpl:end); % Waveform after our signal interval
lennonsigpre  = length(nonsignalpre);
lennonsigpost  = length(nonsignalpost);
noise = nonsignalpre; % Assume that waveform outside our signal interval is pure noise
if lennonsigpost > lennonsigpre
    noise = nonsignalpost; % Choose longer signal
end



%% ************************************************
% Matched filter

%has synch
zt = conv(flipud(pt),yt) * (1/Tsamples); % 1/Tsamples is the dt in the integration (since the integral has been replaced by a sum here)



%% ************************************************
% Calculate sampling  needed, since matched filter above has been time-shifted to make it
% causal

filtertimeoffset = (length(pt) -1);


%% ************************************************
% Sample periodically

zk = zt(1+filtertimeoffset:Tsamples:end);
zk = zk(1:lenXkSym);
%zk start at start of sync


%% ************************************************
% Carrier recovery needed perhaps?

%resynctau = conv(rxsigupsampled,flipud(conj(syncsignal)));

%% ************************************************
% Channel estimation of one-tap channel  zk = h0 xk + nk
% Calculate complex channel gain using pilot symbols
% Gain is estimated only once in packet. Perhaps its should be estimated
% more often?

% zkpilot = zk(1:lenpilotsymbols); % Samples corresponding to pilot symbols
% h0 = (transpose(conj(xkpilot)))*zkpilot/(norm(xkpilot)^2); % Channel gain estimate of one-tap channel.
% disp(['Channel gain is ' num2str(abs(h0)) '*exp(j*' num2str(angle(h0)) ')'])


%% ************************************************
% Equalize one-tap channel using one-tap equalizer

zkPkt = zk(lenSync +1 :end);
vkCoded = []; 
noEq = [];
vk = zk(1:lenSync); %debug
for rnd = 1:numPkts
    zkPart = zkPkt(modPktSize*(rnd-1)+1:modPktSize*rnd);
    zkPlt = zkPart(1:lenPltSym);
    
    h0 = (transpose(conj(xkPlt)))*zkPlt/(norm(xkPlt)^2); 
    if(debug)
            disp(['zkP is: [' num2str(zkPlt.') ']'])
            disp(['xkP is: [' num2str(xkPlt.') ']'])
            disp([' h0 is: [' num2str(h0) ']'])
            disp(['vkP is: [' num2str((zkPlt/h0).') ']'])
    end
    vkPart = zkPart/h0;
    vkMsgTemp = vkPart(lenPltSym+1:end);
    vkPilot = vkPart(1:lenPltSym);
    vkCoded = [vkCoded; vkMsgTemp];
    vk = [vk; vkPart];
    if(debug)
        noEq = [noEq; zkPart(lenPltSym+1:end)]; %debug
    end
end



%% ************************************************
% Decode pilot symbols
% BPSK modulation always used for pilot bits
xkpilothat = sqrt(Ex)*(sign(real(vkPilot)) + 1i*sign(imag(vkPilot))); % Minimum distance detector for 4-QAM

%% Estimate SNRequalizer using pilot symbols
sigmantildesqest = norm(vkPilot(1:lenPltSym)-xkPlt)^2/length(xkPlt);
SNRequalizer = Ex/sigmantildesqest;
disp(['Estimated SNR{equalizer} = ' num2str(10*log10(SNRequalizer)) ' dB, using pilot sequence'])

%% ************************************************
% Decode pilot bits
% BPSK modulation always used for pilot bits
pilotbitshat = (1+xkpilothat)/2 > 0.5; % BPSK

%% ************************************************
% demodulate and decode msg
% Decode message bit using soft decoding of repetition code
% Assume BPSK modulation used for coded bits

if(debug)
    noEqCdd = qamdemod(noEq/d, M, 'gray', 'outputtype', 'bit'); %debug
end

vkCdd = qamdemod(vkCoded/d, M, 'gray', 'outputtype', 'bit');


vkDecoded = [];
cBitsInPkt = 2*msgInPkt;
for rnd = 1:numPkts
    vkBitsPart = vkCdd(cBitsInPkt*(rnd-1)+1:cBitsInPkt*rnd);
    vkDecoded = [vkDecoded; decode(vkBitsPart).'];
end
messagebitshat = vkDecoded;
bitshat = [pilotbitshat; messagebitshat]; % All packet bits




% ************************************************
% Calculate BER of entire packet

BER = mean(messagebitshat(1:length(messagebits)) ~= messagebits);
disp(['BER (excluding pilot) is ' num2str(BER)])
if(debug)
    disp(['MSG is: [' num2str(messagebits(1:msgInPkt).') ']\n']);
    disp(['new is: [' num2str(messagebitshat(1:msgInPkt).') ']\n']);
    disp(['vkD is: [' num2str(vkDecoded(1:msgInPkt).') ']']);
    disp(['rdc is: [' num2str(reDecoded(1:msgInPkt).') ']']);
    disp(['unC is: [' num2str(unCoded(1:msgInPkt).') ']\n']);
    disp('\n--')
    disp('\n--');
    disp(['Cdd is: [' num2str(Coded(1:pktSize).') ']'])
    disp(['vCd is: [' num2str(vkCdd(1:pktSize).') ']'])
    disp(['new is: [' num2str(noEqCdd(1:pktSize).') ']'])
    disp(['xkP is: [' num2str(xkPlt.') ']'])
    disp(['vkP is: [' num2str(vkPilot.') ']'])
    disp(['zkP is: [' num2str(zkPlt.') ']'])
end





% ************************************************
% Plot receiver processing results

thisfignum = 7;

% Timing acquisition using preamble
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(3,1,1)
plot([0:length(syncSignal)-1]/TXusersamplingrate*1e6, real(syncSignal))
ylabel('Ideal')
title('Timing acquisition using preamble')
h(1) = subplot(3,1,2);
plot([0:length(rxsig)-1]/RXusersamplingrate*1e6, real(rxsig))
hold on
plot([0:length(rxsig)-1]/RXusersamplingrate*1e6, imag(rxsig),'r')
ylabel('Received')
axis tight
h(2) = subplot(3,1,3);
plot(([0:length(CTtau)-1]-lenPltSpl+1)/TXusersamplingrate*1e6, abs(CTtau))
ylabel({'abs(C(T,tau))';' '})
axis tight
linkaxes(h,'x')
zoom xon
xlabel('time in microseconds')

% Plot baseband signals (I part)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(2,1,1)
plot([0:length(transmitsignal)-1]/TXusersamplingrate*1e6,real(transmitsignal))
xlabel('time in microseconds')
ylabel('xI(t)')
title('Transmitted baseband signal (real part)')
subplot(2,1,2)
plot([0:length(yt)-1]/(RXusersamplingrate*upsamplingfactor)*1e6,real(yt))
xlabel('time in microseconds')
ylabel('yI(t)')
title('Received and synchronized baseband signal (real part)')



% Plot baseband signals (Q part)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(2,1,1)
plot([0:length(transmitsignal)-1]/TXusersamplingrate*1e6,imag(transmitsignal))
xlabel('time in microseconds')
ylabel('xQ(t)')
title('Transmitted baseband signal (imag part)')
subplot(2,1,2)
plot([0:length(yt)-1]/(RXusersamplingrate*upsamplingfactor)*1e6,imag(yt))
xlabel('time in microseconds')
ylabel('yQ(t)')
title('Received and synchronized baseband signal (imag part)')


% Plot received signal (Frequency response)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
plot(([0:length(transmitsignal)-1]/length(transmitsignal)-0.5)*TXusersamplingrate*1e-6,20*log10(fftshift(abs(fft(transmitsignal)))))
xlabel('frequency in MHz')
ylabel('abs(X(f)) in dB')
title('Received and synchronized baseband signal (frequency response)')


% Plot samples (I part)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(2,1,1)
stem([0:lenXkSym-1],real(xk))
xlabel('time k')
ylabel('xIk')
title('Transmitted symbols (real part)')
subplot(2,1,2)
stem([0:length(vk)-1],real(vk))
xlabel('time k')
ylabel('zIk')
title('Received samples phase corrected (real part)')


% Plot samples (Q part)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(2,1,1)
stem([0:lenXkSym-1],imag(xk))
xlabel('time k')
ylabel('xQk')
title('Transmitted symbols (imag part)')
subplot(2,1,2)
stem([0:length(vk)-1],imag(vk))
xlabel('time k')
ylabel('zQk')
title('Received samples phase corrected (imag part)')


% Signal space (corresponding to pilot symbols)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
plot(real(vkPilot),imag(vkPilot),'o')
hold on
const = sqrt(Ex)*[-1-1i, 1+1i]; 
plot(real(const),imag(const),'rs','MarkerFaceColor','r') % BPSK constellation
xlabel('I')
ylabel('Q')
title('Signal space (Equalized samples corresponding to pilot symbols)')

% Signal space (corresponding to message symbols)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
plot(real(vkCoded),imag(vkCoded),'o')
hold on
plot(real(constellation),imag(constellation),'rs','MarkerFaceColor','r') % BPSK constellation
xlabel('I')
ylabel('Q')
title('Signal space (Equalized samples corresponding to message symbols)')


% Histogram of continuous-time baseband noise
thisfignum = thisfignum+1;
figure(thisfignum)
clf
subplot(2,1,1)
hist(real(noise),50)
xlabel('nI(t)')
title('Histogram of baseband noise n(t)')
subplot(2,1,2)
hist(imag(noise),50)
xlabel('nQ(t)')


% Histogram of discrete-time error (assuming that the equalized channel is  vk = xk + nk)
thisfignum = thisfignum+1;
figure(thisfignum)
clf
nk = vk - xk;
subplot(2,1,1)
hist(real(nk),50)
xlabel('nIk')
title('Histogram of error (assuming that the equalized channel is  vk = xk + nk)')
subplot(2,1,2)
hist(imag(nk),50)
xlabel('nQk')


