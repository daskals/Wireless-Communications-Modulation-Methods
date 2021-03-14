%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Last Version: 23/1/2018
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %The tag has two sensos connected wih a ADC 
  %The packet contains:  
  %Preable              (10 bits)
  %Preable              (2 bits)
  %Sensor ID            (1 bit)
  %Sensor Data section  (10 bits)
  % Each packet section (Preable+Sensor ID+Sensor Data section) is  Hamming
  % Encoded with Hamming(13,19)
  % The tag modoulation is ASK with FM0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [PAM4_enc_packet_out, total_packet_duration, total_packet_length] = Call_Tag_4PAM_Hamming(chanParams)
%% Packets Data 
%  Premble_Data_bits=[1 0 1 0 1 0 1 1 1 1];
%  Tag_ID_bits= randi([0 1],1, 2);
%  Sensor_ID_bits= randi([0 1],1, 1);
%  Sensor_Data_bits= randi([0 1],1, 10);
%  Dummy_bit=1;
 

%% ---------(-3 || 00)----(-1 || 01)------(1 || 11)-----(3 || 10)----------
%% Fixed Data
%(-1, -1, 3,-3,-1)-> 01  01 10 00 01

Premble_Data_bits=[0 1  1 1  1 0  0 0  0 1];
%Premble_Data_bits=[0 1  1 1  1 0  0 0  0 1  1 1   1 0  0 0];
%Premble_Data_bits=[0 0  1 0  1 0  1 0  0 0 ];
Premble_Data_bits=[1 0  0 0  1 0  0 0  1 0  0 1  1 1];

Tag_ID_bits=[0 0];
Sensor_ID_bits=[0 1];
Sensor_Data_bits=[0 1  1 1  1 0  0 0  1 1];

if chanParams.HammingCode==0
    total_packet=[Premble_Data_bits Tag_ID_bits Sensor_ID_bits Sensor_Data_bits];
elseif chanParams.HammingCode==1
    codedIDsDATA=HAmmingcoded18bits([Tag_ID_bits Sensor_ID_bits Sensor_Data_bits]);
    total_packet=[Premble_Data_bits codedIDsDATA];
end
FM0_enc_packet=[];

%% ---------(-3 || 00)----(-1 || 01)------(1 || 11)-----(3 || 10)----------
    j=1;
for i=1:2:length(total_packet)-1
     if total_packet(i)==0 && total_packet(i+1)==0
        %G1 
       FM0_enc_packet(j)=-3;
    elseif total_packet(i)==0 && total_packet(i+1)==1
        %G2
        FM0_enc_packet(j)=-1;
        symb(j)=-1;
    elseif total_packet(i)==1 && total_packet(i+1)==1
        %G3
        FM0_enc_packet(j)=1;
    elseif total_packet(i)==1 && total_packet(i+1)==0 
        %G4
        FM0_enc_packet(j)=+3;
     end
     j=j+1;
end
 
 over=round(chanParams.Tsymbol/chanParams.Ts);
 total_packet_length=length(FM0_enc_packet);
 total_packet_duration=length(FM0_enc_packet)*chanParams.Tsymbol;
 PAM4_enc_packet_out=conv(upsample(FM0_enc_packet, over), ones(over, 1));
 PAM4_enc_packet_out = PAM4_enc_packet_out(1:length(FM0_enc_packet)*over); % truncate
 
end

