% Baseband I,Q modulation to transmit into USRP transmitter
% Mid-project solution
% Only transmits one bit using a repetition code
% Hints for the final part of the project are provided as ** hints **
% As good coding practice, define functions (e.g., bits2symbolMap(), pulseShaping(), etc.) to organize the code logically


% *********************************************************
% Set parameters
clear
clc
set(0,'DefaultTextInterpreter','latex')
disp(' ')

% Fixed (given) parameters. THESE MUST NOT BE CHANGED
filename = 'transmitsignal.mat'; % filename of the USRP transmit signal that will be created.
TXusersamplingrate = 100e6; % user TX baseband sampling rate in Hz. 
maxbandwidth = 25e6; % Maximum user TX RF signal one-sided bandwidth in Hz.
maxlentxsignal = 800e-6; % Maximum length of packet allowed in seconds. 
amplitudebound = 1.0; % Max amplitude allowed. 

% Communication system design parameters
symbolrate = 25e6; % Symbol rate in Hz. Chosen to be a divisor of TXusersamplingrate for ease of implementation. ** Try to increase this substantially for the final project **


debug = false;
b = 2; %doing 2^b-QAM 
M = 2^b;


%for 4-QAM
sigma = 0;
Ex = .0816;
d = sqrt(Ex*((M^2 -1)/12));
dSync = sqrt(2)*d;

% Ex = (d^2/12)(M^2 -1)
% BER = 2(1- 1/M) * Q(sqrt((2*SNR)/(M^2 -1))
% SNR ~ Ex / variance of noise
%snr = 1


lenMsg = 13000; % Total number of message bits. ** Need more bits in final project **
lenPlt = 16; %must be even (len % b = 0)
msgInPkt = 20;
if(mod(lenMsg, msgInPkt) == 0)
else
    disp(['Bad msgInPkt value'])
    return
end
pktSize = msgInPkt + lenPlt;
modPktSize = ((2*pktSize) - lenPlt)/b;
numPkts = lenMsg/msgInPkt;


 % Number of pilot bits (includes timing and frame sync bits). ** Change this for the final project **
lenBuf = 150;
lenSync = 800;

ptRect = 1; 
numinfobits = 2; % k=1 for repetition code. ** Choose 1 if convolutional code used in the final project **
numcodedbits = 1; % n of repetition code. So code rate is k/n. ** Optionally implement a convolutional code in the final project **
alpha = .25;

%% *********************************************************
% Calculate other parameters

maxlentxsignalsamples = ceil(maxlentxsignal*TXusersamplingrate); % Maximum length of packet allowed in samples.
T = 1/symbolrate; % Symbol period in microsec.
Tsamples = ceil(T*TXusersamplingrate); % Symbol period in samples (at user TX baseband sampling rate). The ceil() should not do anything due to choice of T
fs = TXusersamplingrate/symbolrate;

%% *********************************************************
% Choose pulse

% Use rectangular pulse. ** Try to use a bandwidth-efficient pulse for the final project **
if(ptRect)
    pt = ones(Tsamples,1);
    pt = Tsamples*pt/norm(pt); % Nyquist pulse of unit energy. 1/Tsamples is the 'delta t' in the integration
    %pt = pt / (max(abs(pt)));
else
    % Use sqrt-raised cosine filter form  ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as another possible filter
    beta= .1;%Ex;
    span = Tsamples; %nump samples
    samplesPerSymbol  = fs;
    shape = 'sqrt';
    
    pt = rcosdesign(beta, span, samplesPerSymbol, shape); % Tsamples,1/2/T,alpha,fs/T,'rolloff','sqrt');
    
    txfilter = comm.RaisedCosineTransmitFilter;
    txfilter.FilterSpanInSymbols = Tsamples;
    bZ = coeffs(txfilter);
    txfilter.Gain = 1/sum(bZ.Numerator);
    %fvtool(txfilter)
    %return
    %pt = Tsamples*pt/norm(pt);
    pt = txfilter.impz;
    pt = Tsamples*pt/norm(pt);
end

%% *********************************************************
% Constellations

syncConstellation = qammod([0:1], 2, 'gray');
%pilotconstellation = pilotconstellation / (10*max(abs(pilotconstellation)));

pilotConstellation = d*qammod([0:M-1], M, 'gray');
%syncConstellation = syncConstellation / (10*max(abs(constellation)));

constellation = d*qammod([0:M-1], M, 'gray');
%constellation = constellation / (10*max(abs(constellation)));


%r = numinfobits/numcodedbits;
%disp(['Using repetition-coded 4-QAM modulation of code rate r = ' num2str(r) ''])
%disp(['So, bit rate R = ' num2str(r*b) ' (info) bits/symbol'])




%% **********************************************************
% Create packet bits = (pilots, )

% Create message bits
messagebits = (randn(lenMsg,1)) > 0; %actual message
pilotbits = randi(2, lenPlt, 1) -1;  %used for Eq
sbA =(2*randi(2, lenSync/5,1)-1)-2;  %.2
sbB =ones( lenSync/10,1);            %.1
sbC =-1*ones( lenSync/10,1);         %.1
sbD =zeros( lenSync/10,1);           %.1
syncbits = [sbB; sbB; sbA; sbC; flip(sbA); sbC; sbC; sbB];    %timing sync
buffer = zeros(lenBuf, 1);        %so bits aren't lost



%%make packets/frames
allPkts = [];
unCoded = [];
reDecoded = [];
Coded = [];
for rnd = 1:numPkts
    msgPart = messagebits(msgInPkt*(rnd-1)+1:msgInPkt*rnd);
    encodedPart = encode(msgPart).';                     %convolution code
    if(debug)
        unCoded = [unCoded; msgPart]; %debug
        Coded = [Coded; encodedPart]; %debug
        reDecoded = [reDecoded; decode(encodedPart).']; %debug
    end
    allPkts = [allPkts; pilotbits; encodedPart];          %make all pkts
end

%% **********************************************************
% Bits-to-symbol mapping
% Create packet symbols = (pilot symbols, message symbols)

xkPkts = d * qammod(allPkts, M, 'gray', 'inputType', 'bit');
xkPlt  = xkPkts(1:lenPlt/b); %qammod(pilotbits, 4, 'gray', 'inputType', 'bit');
xkSync = dSync * syncbits; %qammod(syncbits, 2, 'gray', 'inputType', 'bit');

lenPktsSym = length(xkPkts);
lenPltSym = length(xkPlt);
lenMsgSym = lenPktsSym - lenPltSym;
lenSyncSym = length(xkSync);

%% **********************************************************
% build xk

xkRel = [xkSync; xkPkts]; %without buffer
xkTot = [buffer; xkRel];

lenXkRel = length(xkRel); 
lenXkTot = length(xkTot); 


%% **********************************************************
% Pulse modulation
% Create I,Q baseband signals (I is real and Q is imaginary component)

% Create the baseband signals by using pulse
% If the signal does not satisfy the allowed maxbandwidth, the transmit
%   mask of the back-end software will filter out the excess bandwidth,
%   thus distorting the signal
xkproposedup = upsample(xkTot,Tsamples); % If Tsamples is not an integer, this method of upsampling will not work
xtproposed = conv(xkproposedup,pt);
lenTotSpl = length(xtproposed); % Length of packet in samples
lenpacket = lenTotSpl/TXusersamplingrate'; % Length of packet in seconds


%% **********************************************************
% Check created transmit signal for length and amplitude limits

flag = true; % Flag to indicate correctness of created transmit signal

% Check transmit signal length falls within allowed limit
if lenpacket > maxlentxsignal
    disp(['Error: Signal not created! Packet length ' num2str(lenpacket*1e3) ' milliseconds should be no more than ' num2str(maxlentxsignal*1e3) ' milliseconds.'])
    beep;
    flag = false;
end
disp(['has size ' num2str(lenpacket*1e3) ' milliseconds'])

% Check transmit signal amplitude falls within allowed limit
mxamplitude = max(max(abs(real(xtproposed))),max(abs(imag(xtproposed))));
if mxamplitude > amplitudebound
    disp(['Error: Signal not created! I and Q amplitude must be no more than ' num2str(amplitudebound)])
    beep;
    flag = false;
    disp(['Found a size of ' num2str(mxamplitude)])
end
disp(['Found a size of ' num2str(mxamplitude)])

%% quit
if ~flag
    % Something is wrong
    return;
end


%% **********************************************************
% Finalize transmit signal and plot it





%% *********************************************************
% Create the USRP transmit signal
% transmitsignal = xI(t)+jxQ(t) sampled at stated sampling rate

transmitsignal = reshape(xtproposed,[],1); % Transmit signal must be column vector

%(transmitsignal.');


%% *********************************************************
% Variables to be used in receiver processing and performance calculation

xkRel = reshape(xkRel,[],1);
xkTot = reshape(xkTot,[],1); % Complex symbols transmitted (sequence may be truncated if too long)
xkPlt = reshape(xkPlt,[],1); % Pilot symbols
xkSync = reshape(xkSync,[],1); % Pilot symbols
pt = reshape(pt,[],1); % Pulse

lenBufSpl = lenBuf*Tsamples; %OUR OWN ZERO BUFFER
lenPltSpl = lenPltSym*Tsamples; % Length of pilot in samples
lenSyncSpl = lenSyncSym*Tsamples;
lenPtSpl = length(pt); % Pulse length in samples
lenTXSpl = length(transmitsignal);

bufSignal = transmitsignal(1:lenBufSpl);%OUR OWN ZERO BUFFER
syncSignal = transmitsignal(lenBufSpl+1: lenBufSpl+lenSyncSpl);
pilotsignal = transmitsignal(lenBufSpl+lenSyncSpl+1:lenBufSpl+lenSyncSpl+lenPltSpl); 

% Calculate number of packet symbols transmitted 
%(including pilot, not buffer though)
lenXkSym = floor((lenTXSpl-lenPtSpl-lenBufSpl+1)/Tsamples); 
% Extract transmitted symbols for performance calculation
xk = xkTot(1:lenXkSym); 



disp(['Symbol period = ' num2str(T*1e6) ' microseconds, which is ' num2str(Tsamples) ' samples'])
disp(['length of packet to be transmitted = ' num2str(lenpacket*1e6) ' microseconds'])
disp(' ')



%% ************************************************
% Save transmit signal, to be uploaded to the USRP
save(filename,'transmitsignal') 



%% ************************************************
% Plots

close all

% Image to be transmitted
if exist('transmitimage','var')
    figure(3)
    imshow(transmitimage,'InitialMagnification','fit')
end


%% Plot time domain signals
ax = [];
figure(4)
clf
subplot(5,1,1);
plot([1:lenPtSpl]/TXusersamplingrate*1e6,pt)
ylabel('p(t)')
xlabel('time  in  microseconds')
title('Time domain signals')
ax1(1) = subplot(5,1,2);
stem([1:lenXkSym],real(xk),'bx')
ylabel('xIk')
ax1(2) = subplot(5,1,3);
stem([1:lenXkSym],imag(xk),'bx')
ylabel('xQk')
xlabel('discrete time  k')
ax2(1) = subplot(5,1,4);
plot([1:lenTotSpl]/TXusersamplingrate*1e6,real(transmitsignal))
ylabel('xI(t)')
ax2(2) = subplot(5,1,5);
plot([1:lenTotSpl]/TXusersamplingrate*1e6,imag(transmitsignal))
ylabel('xQ(t)')
xlabel('time  in  microseconds')
linkaxes(ax1,'x')
linkaxes(ax2,'x')
zoom xon

%% Plot frequency domain signals
figure(5)
clf
subplot(2,1,1)
plot([-lenTXSpl/2+1:lenTXSpl/2]/lenTXSpl*TXusersamplingrate*1e-6,20*log10(abs(fftshift(1/sqrt(lenTXSpl)*fft([pt ; zeros(lenTXSpl-lenPtSpl,1)])))))
ylabel('abs(P(f)) in dB')
vax = axis;
axis([-TXusersamplingrate*1e-6/2 TXusersamplingrate*1e-6/2 vax(3) vax(4)])


xlabel('frequency  in  MHz')
title('Frequency responses')
subplot(2,1,2)
plot([-lenTXSpl/2+1:lenTXSpl/2]/lenTXSpl*TXusersamplingrate*1e-6,20*log10(abs(fftshift(1/sqrt(lenTXSpl)*fft(transmitsignal)))))
ylabel('abs(X{base}(f)) in dB')
vax = axis;
axis([-TXusersamplingrate*1e-6/2 TXusersamplingrate*1e-6/2 vax(3) vax(4)])
xlabel('frequency  in  MHz')



%% Constellation of pilot symbols

maxconstellation = max(abs(constellation));
figure(6)
plot(real(pilotConstellation),imag(pilotConstellation),'rs','MarkerFaceColor','r') % BPSK constellation
hold on;
plot(real(xkPlt),imag(xkPlt),'bo') % BPSK constellation

axis([-1.5*maxconstellation 1.5*maxconstellation -1.5*maxconstellation 1.5*maxconstellation])
xlabel('I')
ylabel('Q')
title('Constellation of pilot symbols')


%% Constellation of message symbols
figure(7)
plot(real(constellation),imag(constellation),'rs','MarkerFaceColor','r')
axis([-1.5*maxconstellation 1.5*maxconstellation -1.5*maxconstellation 1.5*maxconstellation])
xlabel('I')
ylabel('Q')
title('Constellation of message symbols')


%% Baseband I,Q signals created
figure(1)
clf
ax3(1) = subplot(2,1,1);
plot([0:length(transmitsignal)-1]/TXusersamplingrate*1e6,real(transmitsignal))
xlabel('time in microseconds')
ylabel('xI(t)')
title('Transmitted baseband signal')
ax3(2) = subplot(2,1,2);
plot([0:length(transmitsignal)-1]/TXusersamplingrate*1e6,imag(transmitsignal))
xlabel('time in microseconds')
ylabel('xQ(t)')
linkaxes(ax3,'x')
zoom xon


%% Frequency response of complex-valued baseband signal created
figure(2)
clf
plot(([0:length(transmitsignal)-1]/length(transmitsignal)-0.5)*TXusersamplingrate*1e-6,20*log10(abs(fftshift(fft(transmitsignal)))))
xlabel('frequency in MHz')
ylabel('abs(X{base}(f)) in dB')
title('Spectrum of transmitted baseband signal')




    



