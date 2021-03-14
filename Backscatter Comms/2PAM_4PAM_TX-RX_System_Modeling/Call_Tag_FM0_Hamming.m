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
  
function [FM0_enc_packet_out, total_packet_duration, total_packet_length] = Call_Tag_FM0_Hamming(chanParams)

%% Packets Data 
%  Premble_Data_bits=[1 0 1 0 1 0 1 1 1 1];
%  Tag_ID_bits= randi([0 1],1, 2);
%  Sensor_ID_bits= randi([0 1],1, 1);
%  Sensor_Data_bits= randi([0 1],1, 10);
%  Dummy_bit=1;
 

%% Fixed Data
Premble_Data_bits=[1 0 1 0 1 0 1 1 1 1];
Tag_ID_bits=[0 1];
Sensor_ID_bits=[0];
Sensor_Data_bits=[0 1 1 1 1 0 0 0 1 1];
Dummy_bit=1;

if chanParams.HammingCode==0
    total_packet=[Premble_Data_bits Tag_ID_bits Sensor_ID_bits Sensor_Data_bits Dummy_bit];
elseif chanParams.HammingCode==1
    codedIDsDATA=HAmmingcoded18bits([Tag_ID_bits Sensor_ID_bits Sensor_Data_bits]);
    total_packet=[Premble_Data_bits codedIDsDATA Dummy_bit];
end

start=1;
FM0_enc_packet=[];

  for index=1: length(total_packet)

    if total_packet(index)==0 && start==1
        FM0_enc_packet(start)=1;
        FM0_enc_packet(start+1)=0;
     
    elseif total_packet(index)==1 && start==1
        FM0_enc_packet(start)=1;
        FM0_enc_packet(start+1)=1;
        
    elseif total_packet(index)==0 && FM0_enc_packet(start-1)==0
        FM0_enc_packet(start)=1;
        FM0_enc_packet(start+1)=0;
    elseif total_packet(index)==0 && FM0_enc_packet(start-1)==1
        FM0_enc_packet(start)=0;
        FM0_enc_packet(start+1)=1;
        
   elseif (total_packet(index)==1) && FM0_enc_packet(start-1)==1
        FM0_enc_packet(start)=0;
        FM0_enc_packet(start+1)=0;
    elseif (total_packet(index)==1) && FM0_enc_packet(start-1)==0
        FM0_enc_packet(start)=1;
        FM0_enc_packet(start+1)=1;
    end 
    
    start=start+2;

  end
 over=round(chanParams.Tsymbol/chanParams.Ts);
 total_packet_length=length(FM0_enc_packet)/2;
 total_packet_duration=length(FM0_enc_packet)*chanParams.Tsymbol;
 FM0_enc_packet_out=conv(upsample(FM0_enc_packet, over), ones(over, 1));
 FM0_enc_packet_out = FM0_enc_packet_out(1:length(FM0_enc_packet)*over); % truncate

end

