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

Ex = 0.01; % Symbol energy. ** Change this for the final project **

lenmessagebits = 9000; % Total number of message bits. ** Need more bits in final project **
lenpilotbits = 100; % Number of pilot bits (includes timing and frame sync bits). ** Change this for the final project **
lensyncbits = 50;

b = 4; %doing 16-QAM % Since BPSK modulation. ** Must be variable in final project **
numinfobits = 1; % k=1 for repetition code. ** Choose 1 if convolutional code used in the final project **
numcodedbits = 2; % n of repetition code. So code rate is k/n. ** Optionally implement a convolutional code in the final project **
alpha = .5;

% *********************************************************
% Calculate other parameters

maxlentxsignalsamples = ceil(maxlentxsignal*TXusersamplingrate); % Maximum length of packet allowed in samples.
T = 1/symbolrate; % Symbol period in microsec.
Tsamples = ceil(T*TXusersamplingrate); % Symbol period in samples (at user TX baseband sampling rate). The ceil() should not do anything due to choice of T
fs = TXusersamplingrate/symbolrate;

% *********************************************************
% Choose pulse

% Use rectangular pulse. ** Try to use a bandwidth-efficient pulse for the final project **
pt = ones(Tsamples,1); pt = Tsamples*pt/norm(pt); % Nyquist pulse of unit energy. 1/Tsamples is the 'delta t' in the integration
% Use sqrt-raised cosine filter form  ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as another possible filter
%pt = firrcos(Tsamples,1/2/T,alpha,fs/T,'rolloff','sqrt'); 
%pt = Tsamples*pt/norm(pt);

% *********************************************************
% Constellations
% Pilot symbols constellation.
pilotconstellation = sqrt(Ex)*[-1; 1]; % BPSK constellation symbols
% Message constellation.
baseCons = [1.5 - 1i*1.5; 0.5 - 1i*1.5; -1.5 - 1i*1.5; -0.5 - 1i*1.5;
            1.5 - 1i*0.5; 0.5 - 1i*0.5; -1.5 - 1i*0.5; -0.5 - 1i*0.5;
            1.5 + 1i*1.5; 0.5 + 1i*1.5; -1.5 + 1i*1.5; -0.5 + 1i*1.5;
            1.5 + 1i*0.5; 0.5 + 1i*0.5; -1.5 + 1i*0.5; -0.5 + 1i*0.5;];
            
constellation = sqrt(Ex)*baseCons;
r = numinfobits/numcodedbits;
disp(['Using repetition-coded 16-QAM modulation of code rate r = ' num2str(r) ''])
disp(['So, bit rate R = ' num2str(r*b) ' (info) bits/symbol'])




% **********************************************************
% Create packet bits = (sync bits, pilot bits, message bits)

% Create message bits
messagebits = (randn(lenmessagebits,1)) > 0.5;
%messagebits = zeros(lenmessagebits, 1);

% Create pilot bits for channel estimation
% These will also serve for timing recovery 
% These will also serve for frame synchronization because they have a sharp auto-correlation function
pilotbits = (randn(lenpilotbits,1)) > 0.5;
syncbits = (randn(lensyncbits,1)) > 0.5; 
%pilotbits = ones(lenpilotbits, 1);

% Use convolution code
codedbits = encode(messagebits)';

% All bits in frame
bits = [pilotbits; codedbits];



% **********************************************************
% Bits-to-symbol mapping
% Create packet symbols = (pilot symbols, message symbols)

% BPSK modulation for pilot bits. ONLY NEED TO USE BPSK MODULATION FOR PILOT BITS
xkpilot = pilotconstellation(pilotbits+1); % BPSK modulation
%lenpilotsymbols = b*lenpilotbits; 
lenpilotsymbols = 1*lenpilotbits; % Due to BPSK modulation 
xksync = pilotconstellation(syncbits + 1);
%lensyncsymbols = b*lensyncbits;
lensyncsymbols = 1*lensyncbits; % Due to BPSK modulation 

% BPSK modulation for coded bits. ** Use M-QAM in final project **
lencodedsymbols = length(codedbits)/4; %with 16-QAM
toCode = zeros(lencodedsymbols, 1);  %
for i = 1:lencodedsymbols
    one   = codedbits(4*i - 3);
    two   = codedbits(4*i - 2);
    four  = codedbits(4*i - 1);;
    eight = codedbits(4*i - 0);
    toCode(i) = 1*one + 2*two + 4*four + 8*eight;
end

xkmessage = constellation(toCode+1); % BPSK modulation
%xkmessage = OFDMMod(128, 32, xkmessage')';

% Packet symbols = (pilot symbols, message symbols)
%xkproposed = [reshape(xkpilot,[],1); reshape(xkpilot,[],1); reshape(xkmessage,[],1)];
%xkproposed = [reshape(xkpilot,[],1); reshape(xkmessage,[],1)];
%xkproposed = reshape(xkpilot,[],1); % One extra for real transmission
xkproposed = reshape(xkpilot,[],1);
for ii=1:length(xkmessage)/100
    xkproposed = [xkproposed; reshape(xkmessage((ii-1)*100+1:ii*100),[],1)];
        %reshape(xksync,[],1);];
end
xkproposed = [xkproposed xkmessage((ii)*100+1:end)];
lenpacketsymbols = length(xkproposed); % Total number of symbols in packet



% **********************************************************
% Pulse modulation
% Create I,Q baseband signals (I is real and Q is imaginary component)

% Create the baseband signals by using pulse
% If the signal does not satisfy the allowed maxbandwidth, the transmit
%   mask of the back-end software will filter out the excess bandwidth,
%   thus distorting the signal
xkproposedup = upsample(xkproposed,Tsamples); % If Tsamples is not an integer, this method of upsampling will not work
xtproposed = conv(xkproposedup,pt);
lenpacketsamples = length(xtproposed); % Length of packet in samples
lenpacket = lenpacketsamples/TXusersamplingrate'; % Length of packet in seconds


% **********************************************************
% Check created transmit signal for length and amplitude limits

flag = true; % Flag to indicate correctness of created transmit signal

% Check transmit signal length falls within allowed limit
if lenpacket > maxlentxsignal
    disp(['Error: Signal not created! Packet length ' num2str(lenpacket*1e3) ' milliseconds should be no more than ' num2str(maxlentxsignal*1e3) ' milliseconds.'])
    beep;
    flag = false;
end

% Check transmit signal amplitude falls within allowed limit
mxamplitude = max(max(abs(real(xtproposed))),max(abs(imag(xtproposed))));
if mxamplitude > amplitudebound
    disp(['Error: Signal not created! I and Q amplitude must be no more than ' num2str(amplitudebound)])
    beep;
    flag = false;
end



% **********************************************************
% Finalize transmit signal and plot it

if ~flag
    % Something is wrong
    return;
end



% *********************************************************
% Create the USRP transmit signal
% transmitsignal = xI(t)+jxQ(t) sampled at stated sampling rate

transmitsignal = reshape(xtproposed,[],1); % Transmit signal must be column vector



% *********************************************************
% Variables to be used in receiver processing and performance calculation

xkproposed = reshape(xkproposed,[],1); % Complex symbols transmitted (sequence may be truncated if too long)
xkpilot = reshape(xkpilot,[],1); % Pilot symbols
pt = reshape(pt,[],1); % Pulse
lenpilotsamples = lenpilotsymbols*Tsamples; % Length of pilot in samples
lensyncsamples = lensyncsymbols*Tsamples;
lenptsamples = length(pt); % Pulse length in samples
lentxsignalsamples = length(transmitsignal);
pilotsignal = transmitsignal(1:lenpilotsamples); % Pilot signal
syncsignal = transmitsignal(lenpilotsamples+1:lensyncsamples); % sync signal
lenxksymbols = floor((lentxsignalsamples-lenptsamples+1)/Tsamples); % Calculate number of packet symbols transmitted (including pilot)
xk = xkproposed(1:lenxksymbols); % Extract transmitted symbols for performance calculation



disp(['Symbol period = ' num2str(T*1e6) ' microseconds, which is ' num2str(Tsamples) ' samples'])
disp(['length of packet to be transmitted = ' num2str(lenpacket*1e6) ' microseconds'])
disp(' ')



% ************************************************
% Save transmit signal, to be uploaded to the USRP
save(filename,'transmitsignal') 


% ************************************************
% Plots

close all

% Image to be transmitted
if exist('transmitimage','var')
    figure(3)
    imshow(transmitimage,'InitialMagnification','fit')
end


% Plot time domain signals
ax = [];
figure(4)
clf
subplot(5,1,1);
plot([1:lenptsamples]/TXusersamplingrate*1e6,pt)
ylabel('p(t)')
xlabel('time  in  microseconds')
title('Time domain signals')
ax1(1) = subplot(5,1,2);
stem([1:lenpacketsymbols],real(xkproposed),'bx')
ylabel('xIk')
ax1(2) = subplot(5,1,3);
stem([1:lenpacketsymbols],imag(xkproposed),'bx')
ylabel('xQk')
xlabel('discrete time  k')
ax2(1) = subplot(5,1,4);
plot([1:lenpacketsamples]/TXusersamplingrate*1e6,real(transmitsignal))
ylabel('xI(t)')
ax2(2) = subplot(5,1,5);
plot([1:lenpacketsamples]/TXusersamplingrate*1e6,imag(transmitsignal))
ylabel('xQ(t)')
xlabel('time  in  microseconds')
linkaxes(ax1,'x')
linkaxes(ax2,'x')
zoom xon


% Plot frequency domain signals
figure(5)
clf
subplot(2,1,1)
plot([-lenpacketsamples/2+1:lenpacketsamples/2]/lenpacketsamples*TXusersamplingrate*1e-6,20*log10(abs(fftshift(1/sqrt(lenpacketsamples)*fft([pt ; zeros(lenpacketsamples-lenptsamples,1)])))))
ylabel('abs(P(f)) in dB')
vax = axis;
axis([-TXusersamplingrate*1e-6/2 TXusersamplingrate*1e-6/2 vax(3) vax(4)])
xlabel('frequency  in  MHz')
title('Frequency responses')
subplot(2,1,2)
plot([-lenpacketsamples/2+1:lenpacketsamples/2]/lenpacketsamples*TXusersamplingrate*1e-6,20*log10(abs(fftshift(1/sqrt(lenpacketsamples)*fft(transmitsignal)))))
ylabel('abs(X{base}(f)) in dB')
vax = axis;
axis([-TXusersamplingrate*1e-6/2 TXusersamplingrate*1e-6/2 vax(3) vax(4)])
xlabel('frequency  in  MHz')


% Constellation of pilot symbols
figure(6)
plot(sqrt(Ex)*[-1, 1],0,'rs','MarkerFaceColor','r') % BPSK constellation
xlabel('I')
ylabel('Q')
title('Constellation of pilot symbols')


% Constellation of message symbols
maxconstellation = max(abs(constellation));
figure(7)
plot(real(constellation),imag(constellation),'rs','MarkerFaceColor','r')
axis([-1.5*maxconstellation 1.5*maxconstellation -1.5*maxconstellation 1.5*maxconstellation])
xlabel('I')
ylabel('Q')
title('Constellation of message symbols')


% Baseband I,Q signals created
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

% Frequency response of complex-valued baseband signal created
figure(2)
clf
plot(([0:length(transmitsignal)-1]/length(transmitsignal)-0.5)*TXusersamplingrate*1e-6,20*log10(abs(fftshift(fft(transmitsignal)))))
xlabel('frequency in MHz')
ylabel('abs(X{base}(f)) in dB')
title('Spectrum of transmitted baseband signal')




    



