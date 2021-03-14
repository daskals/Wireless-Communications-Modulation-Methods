clc;
close all;
clear all;
%QPSK Transmitter and I-Q Correlation Receiver 
%JC 12/23/05 
%Run from editor debug(F5) 
%m-file for simulating a QPSK transmitter and receiver by modulating with a pseudo 
%random bit stream. A serial to parallel conversion of the pseudo random 
%bit stream is performed with mapping of two bits per symbol(phase). A cosine and 
%sine carrier is configured and the I and Q symbols modulate these 
%carriers. The I and Q carriers are combined and time and frequency domain 
%plots are provided showing key waveforms at various positions in the QPSK 
%transmitter and correlation receiver. A parallel to serial conversion is used on the 
%output of the receiver. The simulation uses a serial "passband" approach. 
%In  other words a carrier is used. Notes and a reference is provided at the end of the m-file. 
%=================================================================== 
clear; 
fcarr=10e3;         % Carrier frequency(Hz) 
N=20;		        % Number of data bits(bit rate) 
fs=40*1e3;		    % Sampling frequency 
Fn=fs/2;            % Nyquist frequency 
Ts=1/fs;	        % Sampling time = 1/fs 
T=1/N;		        % Bit time 
randn('state',0);   % Keeps PRBS from changing on reruns 
td=[0:Ts:(N*T)-Ts]';% Time vector(data)(transpose) 
%=================================================================== 
% The Transmitter. 
%=================================================================== 
data=sign(randn(N,1))';%transpose 
data1=ones(T/Ts,1)*data; 
data2=data1(:); 
 
%display input data bits in command window 
data_2=data2'; 
data_2=data_2 >0; 
x=0; 
transmitted_data_bits=data_2(1:(fs+x)/N:end) 
 
%Serial to parallel (alternating) 
tiq = [0:Ts*2:(N*T)-Ts]';% Time vector for I and Q symbols(transpose) 
 
bs1=data(1:2:length(data));%odd 
symbols=ones(T/Ts,1)*bs1; 
Isymbols=symbols(:);%I_waveform 
 
bs2=data(2:2:length(data));%even 
symbols1=ones(T/Ts,1)*bs2; 
Qsymbols=symbols1(:);%Q_waveform 
 
 
%generate carrier waves 
%cosine and sine wave 
%2 pi fc t is written as below 
twopi_fc_t=(1:fs/2)*2*pi*fcarr/fs;  
a=1; 
%phi=45*pi/180 
phi=0;%phase error 
cs_t = a * cos(twopi_fc_t + phi); 
sn_t = a * sin(twopi_fc_t + phi); 
 
cs_t=cs_t';%transpose 
sn_t=sn_t';%transpose 
si=cs_t.*Isymbols; 
sq=sn_t.*Qsymbols; 
sumiq=si+sq; 
sumiq=.7*sumiq;%reduce gain to keep output at +/- one 
 
%============================================================= 
%Noise 
var=0;%make .1 to 1 to increase noise 
sumiq=sumiq+sqrt(var)*randn(size(sumiq)); 
%============================================================ 
 
%============================================================= 
%Receiver 
%============================================================= 
sig_rx1=sumiq.*cs_t;%cosine 
%simple low pass filter 
rc1=.005;%time constant 
ht1=(1/rc1).*exp(-tiq/rc1);%impulse response 
ycfo1=filter(sig_rx1,1,ht1)/fs; 
 
sig_rx=sumiq.*sn_t;%sine 
%simple low pass filter 
rc=.005;%time constant 
ht=(1/rc).*exp(-tiq/rc);%impulse response 
ycfo=filter(sig_rx,1,ht)/fs; 
 
 
 
bit1=sign(ycfo1);%+/-1 
bit2=sign(ycfo);%+/-1 
bit3=bit1 >0;%0 and 1 
bit4=bit2 >0;%0 and 1 
bitout=[bit3]';%transpose 
bitout1=[bit4]';%transpose 
 
%Parallel to serial bitstream(uses concatenation{joining} and interleaving) 
     
bitout2=[bitout]; 
x=1380;%This is a cluge way to program but x is required to make the parallel 
%to serial converter work if one changes the basic parameters such as N,fs,etc. 
%x=N*(bit3 # 1's or 0's in first bit time)-fs:x=(20*2069)-40000=1380 
bitout2=bitout2(1:(fs+x)/N:end); 
bitout2=[bitout2]; 
bitout3=[bitout1]; 
bitout3=bitout3(1:(fs+x)/N:end); 
bitout3=[bitout3]; 
bitfinalout=[bitout2;bitout3]; 
bitfinalout=bitfinalout(1:end); 
 
%display received output data bits in command window 
Received_data_bits=bitfinalout 
 
%Received data output 
data1a=ones(T/Ts,1)*bitfinalout; 
bitfinal1=data1a(:); 
 
 
 
%===================================================================== 
%Plots 
%====================================================================== 
figure(1) 
subplot(3,2,1) 
plot(td,data2) 
axis([0 1 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('Input Data') 
 
subplot(3,2,3) 
plot(tiq,Isymbols) 
axis([0 1 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('I Channel(one bit/symbol(phase)) Data') 
 
subplot(3,2,5) 
plot(tiq,Qsymbols) 
axis([0 1 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('Q Channel(one bit/symbol(phase)) Data') 
 
subplot(3,2,2) 
plot(tiq,si) 
axis([.498 .502 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('I Channel Modulated Waveform') 
 
subplot(3,2,4) 
plot(tiq,sq) 
axis([.498 .502 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('Q Channel Modulated Waveform') 
 
subplot(3,2,6) 
plot(tiq,sumiq) 
axis([.498 .502 -2 2]); 
grid on 
xlabel('                                                          Time') 
ylabel('Amplitude') 
title('QPSK Output Waveform') 
 
 
%======================================================================== 
%Take FFT of modulated carrier 
%======================================================================== 
y=sumiq; 
NFFY=2.^(ceil(log(length(y))/log(2))); 
FFTY=fft(y,NFFY);%pad with zeros 
NumUniquePts=ceil((NFFY+1)/2);  
FFTY=FFTY(1:NumUniquePts); 
MY=abs(FFTY); 
MY=MY*2; 
MY(1)=MY(1)/2; 
MY(length(MY))=MY(length(MY))/2; 
MY=MY/length(y); 
f1=(0:NumUniquePts-1)*2*Fn/NFFY; 
%========================================================================= 
%Plot frequency domain 
%========================================================================= 
figure(2) 
subplot(3,1,1); plot(f1,MY);xlabel('');ylabel('AMPLITUDE'); 
axis([9500 10500 -.5 1]);%zoom in/out 
title('Frequency domain plots'); 
grid on 
 
subplot(3,1,2); plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(Hz)');ylabel('DB'); 
axis([9000 11000 -80 10]);%zoom in/out 
grid on 
title('Modulated QPSK carrier') 
 
 
figure(3) 
subplot(3,2,1); 
plot(td,bitfinal1) 
title('Received output data'); 
grid on; 
axis([0 1 -2 2]); 
 
subplot(3,2,3); 
plot(tiq,ycfo1); 
title('Filtered I Channel Data'); 
grid on; 
 
subplot(3,2,5); 
plot(tiq,ycfo); 
title('Filtered Q Channel Data'); 
grid on; 
 
subplot(3,2,2); 
plot(tiq,sig_rx1); 
grid on; 
title('Unfiltered I Channel Output'); 
 
subplot(3,2,4); 
plot(tiq,sig_rx); 
grid on; 
title('Unfiltered Q Channel Output'); 
 
phasevl=atan2(ycfo1,ycfo); 
subplot(3,2,6); 
plot(tiq,phasevl); 
grid on; 
title('Output phase voltage levels') 
 
 
%NOTE 
%Serial to parallel conversion of a serial bit stream and mapping of 
%two bits to a symbol(phase) can sometimes be confusing. I will try and explain 
%with an example. 
%Suppose you have a serial bit stream of ten  0 0 1 1 0 1 1 0 1 1 even # of bits       
                                %odd bits     0   1   0   1   1 
                                %even bits      0   1   1   0   1 
                                 
%The odd bits are the I Channel Data at one half the original serial bit stream 
%bit rate. Notice that the amplitudes are +/- one as shown in figure 1. 
%The possible combinations are -1 -1, 1 1, -1 1, 1 -1(four phases or four symbols). 
%The amplitudes, in theory, should be held at +/- 0.707 to keep the summed output of the  
%QPSK transmitter at a  constant amplitude of +/- one. 
%The even bits are the Q Channel Data at one half the original serial bit 
%stream bit rate. Same info as above. 
 
%Things to do: 
%Implement BER code to prove that the BER of the output of either the I or 
%Q channel is equal to the BER of BPSK. Also prove that the symbol 
%BER(combined output) is ~ 3DB poorer than the I or Q channel output. 
%Implement Grey coding and prove that the BER of each approaches equality 
%at high S/N ratios. 
%Implement different types of low pass filters for best BER. Appropriate TX and RX 
%bandpass filters could also be added for best BER. 
 
%A good reference discussing a QPSK Transmitter and look up tables and Gray 
%coding can be found at http://cnx.rice.edu/content/m10042/latest/ 

