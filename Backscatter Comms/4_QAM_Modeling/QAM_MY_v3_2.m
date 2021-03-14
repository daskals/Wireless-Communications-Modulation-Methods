%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis
clc;
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pause(6)
%%RTL SDR parameters
GAIN=-15;
F_ADC = 1e6;
DEC = 1;
Fs = F_ADC/DEC;
Ts = 1/Fs;

%% Sympol parameters
Tsymbol = 1e-3; %ms Rymbol Rate=1/Tymbol Bit Rate= 2/Tymbol
Tbit=Tsymbol;
over = round(Tsymbol/Ts);
newover = 10;



%% PAcket parameters
preamble_length=10;                % NoFM0_prample=[1 0 1 0 1 0 1 1 1 1];
id_length=2;                       % NoFM0_ID=[0 1];
util_length=2;                     % NoFM0_util=[0 1];
codeword_length=12;                % NoFM0_DATA=[0 0 1 1 1 1 0 0 0 1 0 1];
dummybit=0;

total_packet_length=id_length+preamble_length+util_length+codeword_length+dummybit;
total_packet_duration=(total_packet_length)*Tbit;


% Create the preamble in FM0 format.
packet=[1 0 1 0 1 0 1 1 1 1  0 1  0 1  0 0 1 1 1 1 0 0 0 1 0 1 ];


%% Sigmal Prosesing  Variables
Resolution = 1;   % in Hz
N_F = Fs/Resolution;
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;

framelength=3;
t_sampling = framelength*total_packet_duration;     % Sampling time frame (seconds).
N_samples = round(Fs*t_sampling);
t = 0:Ts:t_sampling-Ts;


%%
T = Tsymbol;                              % symbol period
A = 8;                              % filter delay
beta = 0.9;                         % rolloff facto
%%% use MATLAB's built-in function to create a square-root rcosine
g_T = rcosine(1/T, Fs, 'sqrt', beta, A)/sqrt(T);

%QAM seperation

oddbits= packet(1:2:length(packet)); %odd
evenbits=packet(2:2:length(packet)); %even

x_tag_odd= conv(upsample(oddbits, over), ones(over, 1));
x_tag_odd= x_tag_odd(1:length(oddbits)*over);


x_tag_even= conv(upsample(evenbits, over), ones(over, 1));
x_tag_even= x_tag_even(1:length(evenbits)*over);

%N=length(packet);
N=length(oddbits);


%%

%  G1 = unifrnd(0, 0.25)*exp(j*1);
%  G2 = unifrnd(0.26, 0.5)*exp(j*angle(G1));
%  G3 = unifrnd(0.51, 0.75)*exp(j*angle(G1));
%  G4 = unifrnd(0.76, 1)*exp(j*angle(G1));
 
%      % QAM reflection coefficients
 G1 = 1*exp(j*pi/4);
 G2 = 1*exp(j*angle(G1) + j*2*pi/4);
 G3 = 1*exp(j*angle(G1) + j*4*pi/4);
 G4 = 1*exp(j*angle(G1) + j*6*pi/4);

 for ind=1:10

     if  (oddbits(ind)==1 & evenbits(ind)==1)
        %G1 
        preamble(ind)=G1;
    elseif (oddbits(ind)==0 & evenbits(ind)==1)
        %G2
        preamble(ind)=G2;
    elseif (oddbits(ind)==0 & evenbits(ind)==0)
        %G3
        preamble(ind)=G3;
    elseif (oddbits(ind)==1 & evenbits(ind)==0)    
        %G4
        preamble(ind)=G4;
     end     
end

 
 
for i=1:length(x_tag_odd)
 
     if  (x_tag_odd(i)==1 && x_tag_even(i)==1)
        %G1 
        x_G(i)=G1;
    elseif (x_tag_odd(i)==0 && x_tag_even(i)==1)
        %G2
        x_G(i)=G2;
    elseif (x_tag_odd(i)==0 && x_tag_even(i)==0)
        %G3
        x_G(i)=G3;
    elseif (x_tag_odd(i)==1 && x_tag_even(i)==0)    
        %G4
        x_G(i)=G4;
     end
end




%% create training sequence
N_tr = 10;                          % no of training symbols
train_step = 1;                     % spacing of training symbols
train_seq = ones(1, N_tr) + j*ones(1, N_tr);    % training sequence of (1+j)'s
tr_ind = [0:N_tr-1]*train_step + 1;          % training symbols indices
a_n(tr_ind) = train_seq;            % put the trainings in the symbol sequence


train_seq=preamble;

%%


figure (1);
plot(x_G, 'o');
title('TAG symbols');

display_figures=1;

 %% carrier
F_sin = 15e3;            % rate of sine
m = cos(2*pi*F_sin*t);   % baseband sine

F_sin1=19e3;
F_sin2=38e3;
F_sin3=57e3;
m2 = cos(2*pi*F_sin1*t)+cos(2*pi*F_sin2*t)+cos(2*pi*F_sin2*t);   % baseband sine

x_FM = NaN*ones(size(m));   % modulated FM signal
f_c = 0;            % carrier (baseband = 0)
CFO = 10;         % frequency offset
Df = 75e3;          % freq deviation



for tt=1:length(t)
    x_FM(tt) = exp(j*2*pi*(f_c+CFO)*t(tt) + j*2*pi*Df*sum(m2(1:tt))*Ts);
end

    
                    % fft
                   x_fft = fftshift(fft(x_FM, N_F));
                   F_sensor_est_power=10*log10((abs(x_fft).^2)*Ts/50*1e3)-15; 
                   figure(2);
                   subplot(2,1, 1)
                   plot(F_axis/1000000, F_sensor_est_power);
                   title('Frequency Domain of Source signal')  
                   ylabel('Power (dBm)');
                   xlabel('Frequency (MHz)');
                   drawnow;  
                   subplot(2,1, 2)
                   time_axis= 0:Ts:Ts*length(x_FM)-Ts; 
                   plot(time_axis,x_FM);
                   title('Time Domain of Source signal')  
                   xlabel('Time (Sec)');
                   
                 

 %Modulation of signal                                 
L_padding = 20000;
x3_packet= [0*ones(1, L_padding) x_G  0*ones(1, length(x_FM)-length(x_G)-L_padding)];
x_R = x_FM.*x3_packet; % full signal


                   
                   figure(3);
                   subplot(2, 1, 1);
                    plot(time_axis,x_FM);
                    title('Time Domain of Source signal ')  
                    xlabel('Time (Sec)');
                  
                   subplot(2, 1, 2);
                    plot(time_axis, x_R);
                    title('Time Domain of Modulated signal ')  
                    xlabel('Time (Sec)');
                  
            
             scatterplot(x_R)
             title('Symbols Modulated with FM signal ') 
             %%
             % LPF
             
LPF_cutoff=1000;
LPF_upfreq=20000;
LPF_cutfreq=18000;
x_R_fft = fftshift(fft(x_R));           % fft(x_R)
LPF = zeros(size(x_R_fft));             % create ideal LPF
N_F = length(LPF);                      %
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;      %
%LPF(F_axis <= LPF_upfreq & F_axis>= LPF_cutfreq) = 1;     % cutoff freq = 10Hz
LPF(abs(F_axis)<= LPF_cutoff) = 1;     % cutoff freq = 10Hz
x_R_filt_fft = x_R_fft.*LPF;            % filter in freq domain
x_R_filt = ifft(ifftshift(x_R_filt_fft)); % ifft
%x_R_filt=x_R;

                   x_fft = fftshift(fft(x_R_filt, N_F));
                   F_sensor_est_power2=10*log10((abs(x_fft).^2)*Ts/50*1e3)-15; 
                   figure(11);
                   plot(F_axis/1000000, F_sensor_est_power2);
                   title('Frequency Domain  LPF')  
                   ylabel('Power (dBm)');
                   xlabel('Frequency (MHz)');
                   drawnow;  
                   

             scatterplot(x_R_filt)
             title('Symbols Modulated after  LPF ') 

%% coarse CFO estimation/cancellation
N_F = length(x_R_filt);
f_axis = -1/2:1/N_F:(1/2)-(1/N_F);

x_R4_fft = abs(fftshift(fft(x_R_filt.^4, N_F)));

[maxval_x_R4_fft, maxpos_x_R4_fft] = max(x_R4_fft);
coarse_CFO_est = f_axis(maxpos_x_R4_fft)/4;
t_x_discrete = (0:length(x_R_filt)-1);
x = x_R_filt.*exp(-j*2*pi*coarse_CFO_est*t_x_discrete);




          
             scatterplot(x)
             title('Symbols Modulated after CFO'); 
             
        
           %% Filtering
            g_R = g_T;                          % RX filter is matched to TX filter
            x_low = x;
            x_CRRC = conv(x, g_R)*Ts;            % filter x(t) with SRRC

             
           
            scatterplot(x_CRRC)
            title('Symbols Modulated after CRRC filtering'); 
             
  
            %% packet sync using energy
A=8;
d1 = 1;             % first sample to calculate energy
d2 = 3*A*over;      % last sample to calculate

E = zeros(1, d2);    % energies vector

for d = d1:d2
    E(d) = norm(x_low(d:over:d+(N-1)*over))^2;
     %energy_synq(d)=sum(abs(x_low(k : k+total_packet_length*2*newover-1)).^2);
end

if display_figures
    figure(8);
    plot([d1:d2], E(d1:d2));
    title('E(d)');
    xlabel('d (samples)');
    ylabel('EY_l_o_w^2');
    axis tight;
    grid on;
end

% find max energy and d*
[max_x_low pos_max_x_low] = max(E);
pos_max_x_low
%% get samples from x_low_d*(t)
r_k = x_low(pos_max_x_low:over:pos_max_x_low + (N-1)*over);

if display_figures
    
    figure(9);
    subplot(2, 1, 1);
    stem(real(x3_packet), 'r');
    title('In-phase sampled symbols T');
    % plot imaginary part of sampled symbols r_k
    subplot(2, 1, 2);
    stem(real(r_k), 'r');
    title('In-phase sampled symbols R');
    
    % plot real part of sampled symbols r_k
    figure(10);
    subplot(2, 1, 1);
    stem(imag(x3_packet), 'r');
    title('Quadrature sampled symbols T');
    % plot imaginary part of sampled symbols r_k
    subplot(2, 1, 2);
    stem(imag(r_k), 'r');
    title('Quadrature sampled symbols R');
end
      

