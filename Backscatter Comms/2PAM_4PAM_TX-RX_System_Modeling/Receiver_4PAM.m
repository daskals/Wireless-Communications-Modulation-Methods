
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [total_envelope, packets,index, Bit_errors, fixedpacketdata_len, returnthress] = Receiver_4PAM(r, packets, total_packet_length,chanParams)


Thressreceiver=chanParams.Thressreceiver;


Fs = 1e6;

Ts=1/Fs;

newover = 10;
over = round(chanParams.Tsymbol/Ts);

Bit_errors=0;
%% Sigmal Prosesing  Variables
Resolution = 1;   % in Hz
N_F = Fs/Resolution;
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;

DEBUG_en1=0;
DEBUG_en2=1;

% Preamble in FM0 format with symbols (not bits).
preamble=[-1, -1, 3, -3, -1];

preamble_neg=-1*preamble;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decoder  FM0 vectors
bits_FM0_2sd_wayB=[]; 
decision_bits_B=[];
final_packet=[];
returnthress=0;

fixedpacketdata=[0 1  0 1  0 1 1 1 1 0 0 0 1 1];  % id + sensor_id + fixedata  
fixedpacketdata_len=length(fixedpacketdata);
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
           x = r.';  % receive signal
          
           packets = packets + 1
           %fprintf('Packet Window=%d|\n',packets)
           
           %% Absolute operation removes the unknown CFO
            abstream=abs(x).^2;
           %% Matched filtering
            matcheds=ones(round(chanParams.Tsymbol/Ts),1); % the pulse of matched filter has duration Tsymbol
            dataconv=conv(abstream,matcheds);   %  aply the filter with convolution
      
            
           %% Signal Smoothing with  1-D median filtering
           %dataconv= medfilt1(dataconv); 

            %% Downsample same prosedure
            total_env_ds = dataconv(1:over/newover:end); %% by factor of 10 to reduce the computational complexity
            
            %% Time sync of downsample
            total_envelope = total_env_ds(newover+1:end-(newover+1)); % total_env_ds(newover+1:end-newover+1); 
            
            %% remove the DC offset
            total_envelope=total_envelope-mean(total_env_ds);
            
            %scale/normalize values in a vector to be between -3 and 3
            %total_envelope=total_envelope(:)*sqrt(length(total_envelope))/norm(total_envelope);
            total_envelope = -3 + 6.*(total_envelope - min(total_envelope))./(max(total_envelope) - min(total_envelope));
             %% Print Plots

            %averagepower
             averagepower= sum( mean(total_envelope.^2) )/length(total_envelope) 
           if DEBUG_en1==1;
                time_axis= 0:Ts:Ts*length(abstream)-Ts;        %same as xaxis_m= (1: length(abstream))*Ts Captured signal time axis.
                     % fft
                     x_fft = fftshift(fft(x, N_F));
                     F_sensor_est_power=10*log10((abs(x_fft).^2)*Ts/50*1e3)-15; 
                 figure(3);
                  subplot(2, 1, 1);
                    plot(time_axis,abstream);
                    title('Absolute-squared', 'FontSize',14 )  
                    xlabel('Time (Sec)', 'FontSize',12, 'FontWeight','bold');
                    ylabel('Amplitude', 'FontSize',12, 'FontWeight','bold');
                    grid on;
                   subplot(2, 1, 2);
                    plot(F_axis/1000000, F_sensor_est_power);
                    title('Frequency Domain')  
                    xlabel('Frequency (MHz)');
                   drawnow;               
           end
           
            if DEBUG_en2==1;
                 figure(4);
                    subplot(2, 1, 1);
                    plot(dataconv);
                    title('Matched-filtered' ,'FontSize',14 )
                    xlabel('Time (Sec)', 'FontSize',12, 'FontWeight','bold');
                    ylabel('Amplitude', 'FontSize',12, 'FontWeight','bold');
                    grid on;
                    subplot(2, 1, 2);
                    plot(total_envelope);
                    title('DOWNSAMPLED DC Removal')
                    drawnow;                 
            end 

            
              %%  Position estimation of packet with  packet's energy synchronization
%             for k=1:1: length(total_envelope)-(total_packet_length*2*newover)+1
%                     energy_synq(k)=sum(abs(total_envelope(k : k+total_packet_length*2*newover-1)).^2);          
%             end
%             % find the starting point of packet
%             [returnthress  energy_sinq_ind]=max(energy_synq);
%             
%            if returnthress < chanParams.Thressreceiver                          %|| energy <= ENERGYTHRESS/3
%                index=-1;
%                disp '-----Drop packet (NOISE)------';
%              return ;
%            end 
            
           
            %% create the preamble neover format
            preample_neover=upsample(preamble, newover);
            preample_neg_neover=upsample(preamble_neg, newover);
            
            %% Sync via preamble correlation
            corrsync_out = xcorr(preample_neover, total_envelope);
            corrsync_out_neg = xcorr(preample_neg_neover, total_envelope);
            
            [m ind] = max(corrsync_out);
            [m_neg ind_neg] = max(corrsync_out_neg);
             %notice that correlation produces a 1x(2L-1) vector, so index must be shifted.
             %the following operation points to the "start" of the packet.

             
            if (m < m_neg)
               start = length(total_envelope)-ind_neg;
               total_envelope=-total_envelope;
            else
               start = length(total_envelope)-ind;
            end
            
            if(start <= 0)
                index=1;
                disp 'Negative start';
                return ;
            elseif start+((total_packet_length))*newover > length(total_envelope)  %% Check if the detected packet is cut in the middle.
                index=2;
                disp 'Packet cut in the middle!';
                return ;
            end 
            shifted_sync_signal_B=total_envelope(start+length(preample_neover)+1: start+total_packet_length*newover);
            midle_symbol_points=shifted_sync_signal_B(1:newover:end);
          
           
%             
%             figure(6);
%             plot(shifted_sync_signal_B);
%             drawnow;
             figure(7);
             plot(midle_symbol_points,'o');
             drawnow;
%             
       % quantize the input signal x to the alphabet
       % using nearest neighbor method
       averagepower
       alphabet=[]
            alphabet=[-3; -1; 1; 3];
            x=midle_symbol_points(:);
            alpha=alphabet(:,ones(size(x)))'; %matrix 7x4 
            dist=(x(:,ones(size(alphabet)))-alpha).^2; %matrix 7x4 with distances
            [v,iv]=min(dist,[],2); % find the minimum distance in each line .. each column coresponds to a symbol (-3 -1 1 3)
            y_bits=alphabet(iv);    
           
            bitsind=1;
            %(-3 || 00)----(-1 || 01)------(1 || 11)-----(3 || 10)
            
            for i=1:length(y_bits)
                if y_bits(i)==-3
                    decision_bits_B (bitsind)=0;
                    decision_bits_B (bitsind+1)=0;
                elseif y_bits(i)==-1
                    decision_bits_B (bitsind)=0;
                    decision_bits_B (bitsind+1)=1;
                elseif y_bits(i)==1
                    decision_bits_B (bitsind)=1;
                    decision_bits_B (bitsind+1)=1;
                elseif y_bits(i)==3
                    decision_bits_B (bitsind)=1;
                    decision_bits_B (bitsind+1)=0;
                end
                bitsind=bitsind+2;
            end 
                

             if chanParams.HammingCode==1  
              [final_packet] = HAmmingDecoded18bits(decision_bits_B);
            else
                final_packet=decision_bits_B;
             end
             
                if  isequal(final_packet, fixedpacketdata)
                      disp 'Packet Correct !!!!!!!!!!!!!!!!!!!!!!!!!';
                       index=3;        
                else
 
                      disp 'Packet WRONGGGG------------------------';
                      index=4;   
                      Bit_errors = sum(xor(final_packet,fixedpacketdata));
                end  

end

