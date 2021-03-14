%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Last Version: 31/1/2018
% Matlab Vesrion: R2017b (Academic Use) 
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 
close all; 
clear all;

%% Main Program

%% FM Emitter parameters
   chanParams.F_sin = 15e3;               % rate of sine
   chanParams.f_c = 0;                    % carrier (baseband = 0)
   chanParams.CFO = 0;                    % frequency offset
   chanParams.Df = 75e3;                  % freq deviation
   
%% Reader Parameters
    chanParams.Fs = 1e6;
    chanParams.Ts=1/chanParams.Fs;
    chanParams.framelength=3; 
    chanParams.Thressreceiver=200;

%% Tag parameters
    %chanParams.Tsymbol = 5.6e-3; %Tsymbol of the tag
    chanParams.Tsymbol = 5.8e-6; %Tsymbol of the tag
    chanParams.HammingCode=0;   %Hamming code enabled
    
%% Noise Param
   chanParams.noiseSNR=20;
   chanParams.P_TX = 20; 
   
%% Modulation Selection Param
   chanParams.Modtype=0; % 0 for 2PAM ----1 for 4 PAM
   
%% Initiallization
 %Plots
   plot1=1;
   plot2=1;
 % packets status indicator  
   drop_packets=0;
   correct_packets=0;
   error_packets=0;
   negative_starts=0;
   cut_packets=0;
 % other   
   Bit_errors_sum=0;
   packets=0;
   fixedpacketdata_len=0;
 %% PacketsSimulation
    N=100;
  
  
if chanParams.Modtype==0
[FM0_enc_packet, total_packet_duration ,total_packet_length] = Call_Tag_FM0_Hamming(chanParams);
else
[enc_packet, total_packet_duration, total_packet_length] = Call_Tag_4PAM_Hamming(chanParams);
end

[x_FM]= FM_carrier_channel(total_packet_duration, chanParams);
  
 t_sampling = chanParams.framelength*total_packet_duration;     % Sampling time frame (seconds).
 N_samples = round(chanParams.Fs*t_sampling);
 t = 0:chanParams.Ts:t_sampling-chanParams.Ts;

for (steps=1:1:N)
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Tag Part (Transmitter) 
  if chanParams.Modtype==0
  [enc_packet, total_packet_duration, total_packet_length] = Call_Tag_FM0_Hamming(chanParams);
  else
  [enc_packet, total_packet_duration, total_packet_length] = Call_Tag_4PAM_Hamming(chanParams);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if plot1==1
  figure(1);
  plot(enc_packet)
  title('Wavefor for RF transistor')
  xlabel('Time (sec)');
  ylabel('Amplitude (a.u)');
  drawnow;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Channel and Modulation Part
  [r]=modulation_signal(enc_packet, x_FM, chanParams, t);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

  time_axis= 0:chanParams.Ts:chanParams.Ts*length(r)-chanParams.Ts;
  if plot2==1
  figure(2);
  plot(time_axis,abs(r).^2)
  title('FM Modulated signal')
  xlabel('Time (sec)');
  ylabel('Amplitude (a.u)');
  drawnow;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Receiver Part
  if chanParams.Modtype==0
  [total_envelope, packets, index, Bit_errors, fixedpacketdata_len] = Receiver_ASK_FM0(r,packets, total_packet_length, chanParams);
  else
  [total_envelope, packets,index, Bit_errors, fixedpacketdata_len]= Receiver_4PAM_v2(r,packets, total_packet_length, chanParams);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (index==-1) 
      drop_packets=drop_packets+1; 
     %return;
  elseif (index==1) 
      negative_starts=negative_starts+1; 
      %return;
  elseif(index==2)
      cut_packets = cut_packets + 1;
      %return;
  elseif(index==3)
      correct_packets=correct_packets+1;
      %return;
  elseif(index==4)
      error_packets=error_packets+1;
      Bit_errors_sum=Bit_errors_sum+Bit_errors;
      %return;
  end
end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Print Results 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        fprintf('\n');
        fprintf('Emitter: \n');
        fprintf('CFO=%d| Frequency Deviation=%d\n',chanParams.CFO, chanParams.Df)
        fprintf('------------------------------------\n')
        fprintf('Tag: \n');
        fprintf('HammingEncoded=%d|Sended Packets=%d|Sended Bits=%d\n', chanParams.HammingCode, steps, (correct_packets+error_packets)*fixedpacketdata_len)
        fprintf('------------------------------------\n')
        fprintf('Receiver: \n')
        fprintf('Drop Packets=%d\n',drop_packets)
        fprintf('Corecct Packets=%d|Packet Errors=%d\n',correct_packets, error_packets)
        fprintf('Negative Starts=%d|Cut Packets=%d\n', negative_starts, cut_packets)
        %PER is the number of incorrectly received data packets divided by the total number of received packets.
        fprintf('(PER) Packet Error Rate=%d \n', error_packets / (correct_packets+error_packets))
        fprintf('(BER) Bit Error Rate=%d \n', Bit_errors_sum/((correct_packets+error_packets)*fixedpacketdata_len))
  