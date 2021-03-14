%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEL303 - Communication Systems II, Spring 2012       %
% ECE Dept, Technical University of Crete              %
% Project: QPSK Modem                                  %
% John Kimionis (jkimionis@gmail.com)                  %
% Instructor: Aggelos Bletsas (aggelos@telecom.tuc.gr) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;


display_figures = 1;                % turn figures on/off
add_noise = 0;                      % turn noise on/off

Fc = 20;                            % carrier frequency
DF = (12.4/100)*Fc;                 % carrier frequency offset
Dphi = unifrnd(-pi, pi);            % carrier phase offset

LPF_cutoff = 10;                    % RX filter bandwidth
if add_noise
    var_w = 0.0001;                       % noise variance
else
    var_w = 0;                       % noise variance
end
%% transmit pulse parameters
T = 1;                              % symbol period
Ts = 1/100;                         % sampling period
over = T/Ts;                        % oversampling factor
Fs = 1/Ts;                          % sampling frequency
A = 6;                              % filter delay
beta = 0.5;                         % rolloff factor

%%% use MATLAB's built-in function to create a square-root rcosine
g_T = rcosine(1/T, Fs, 'sqrt', beta, A)/sqrt(T);

%% create 4-QAM symbols
N = 100;                             % no of QAM symbols
a_n = sign(randn(1, N)) + j*sign(randn(1, N));
scatterplot(a_n) 
%% create training sequence
N_tr = 20;                          % no of training symbols
train_step = 5;                     % spacing of training symbols
train_seq = ones(1, N_tr) + j*ones(1, N_tr);    % training sequence of (1+j)'s
tr_ind = [0:N_tr-1]*train_step + 1;          % training symbols indices
a_n(tr_ind) = train_seq;            % put the trainings in the symbol sequence

%% create baseband QAM waveform
s_c = conv(upsample(a_n, over), g_T);
t_s = -A*T:Ts:N*T + A*T-Ts;

%% PSD of baseband waveform
N_F = length(s_c);
P_s_c = abs(fftshift(fft(s_c, N_F))*Ts).^2/(Ts*length(s_c));
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;

if display_figures
    figure;
    semilogy(F_axis, P_s_c);
    grid on;
end

%% create passband transmitted waveform
s = real(s_c.*exp(j*2*pi*Fc*t_s));

%% PSD of passband waveform
P_s = abs(fftshift(fft(s, N_F))*Ts).^2/(Ts*length(s));

if display_figures
    figure;
    semilogy(F_axis, P_s);
    grid on;
end

%% channel
h_0 = randn;                            % random channel tap
h = zeros(1, over);                     % create channel with one tap
h(floor(unifrnd(1, over))) = h_0;
t_h = 0:Ts:Ts*length(h)-Ts;             % channel time axis

y = conv(s, h);                       % convolve with channel
t_y = min(t_s)+min(t_h):Ts:max(t_s)+max(t_h);

w = sqrt(var_w)*randn(size(t_y));       % create noise
y = y + w;                              % add noise to channel output

%% demodulate
x_RI = y.*cos(2*pi*(Fc+DF)*t_y + Dphi); % In-phase demodulated component
x_RQ = -y.*sin(2*pi*(Fc+DF)*t_y + Dphi);% Quadrature demodulated component

x_R = x_RI + j*x_RQ;                    % complex demodulated signal

% LPF
x_R_fft = fftshift(fft(x_R));           % fft(x_R)
LPF = zeros(size(x_R_fft));             % create ideal LPF
N_F = length(LPF);                      %
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;      %
LPF(abs(F_axis) <= LPF_cutoff) = 1;     % cutoff freq = 10Hz
x_R_filt_fft = x_R_fft.*LPF;            % filter in freq domain
x_R_filt = ifft(ifftshift(x_R_filt_fft)); % ifft
x_R = x_R_filt;

%% coarse CFO estimation/cancellation
N_F = length(x_R_filt);
f_axis = -1/2:1/N_F:(1/2)-(1/N_F);

x_R4_fft = abs(fftshift(fft(x_R_filt.^4, N_F)));

[maxval_x_R4_fft, maxpos_x_R4_fft] = max(x_R4_fft);
coarse_CFO_est = f_axis(maxpos_x_R4_fft)/4;
t_x_discrete = (0:length(x_R_filt)-1);
x = x_R_filt.*exp(-j*2*pi*coarse_CFO_est*t_x_discrete);

%% Filtering
g_R = g_T;                          % RX filter is matched to TX filter

x_low = conv(x, g_R)*Ts;            % filter x(t) with SRRC

%% packet sync using energy
d1 = 1;             % first sample to calculate energy
d2 = 3*A*over;      % last sample to calculate

E = zeros(1, d2);    % energies vector

for d = d1:d2
    E(d) = norm(x_low(d:over:d+(N-1)*over))^2;
end

if display_figures
    figure;
    plot([d1:d2], E(d1:d2));
    title('E(d)');
    xlabel('d (samples)');
    ylabel('EY_l_o_w^2');
    axis tight;
    grid on;
end

% find max energy and d*
[max_x_low pos_max_x_low] = max(E);

%% get samples from x_low_d*(t)
r_k = x_low(pos_max_x_low:over:pos_max_x_low + (N-1)*over);

if display_figures
    % plot real part of sampled symbols r_k
    figure;
    subplot(2, 1, 1);
    stem(real(r_k), 'r');
    title('In-phase sampled symbols');
    % plot imaginary part of sampled symbols r_k
    subplot(2, 1, 2);
    stem(imag(r_k), 'r');
    title('Quadrature sampled symbols');
end

scatterplot(r_k);
%% freq offset estimation
% create z_l
% r(tr_ind) are the received symbols which correspond to the training
% symbols a_n(tr_ind)
% create z = received_training*conj(training)
z = r_k(tr_ind).*conj(train_seq);   % = r_k(tr_ind).*a_n(tr_ind)';

%% calculate the fourier transform of z and find max
N_F = 1024;             % points to calculate the fourier transform
dF = 1/N_F;             % frequency axis step
f_axis = -1/2:dF:(1/2)-dF;    % frequency axis

Z = fftshift(abs(fft(z, N_F)*Ts));

if display_figures
    % plot the fourier transform of z
    figure;
    plot(f_axis, Z);
    axis tight;
    title('Fourier Transform of z_l');
    xlabel('f (Hz)');
    ylabel('|Z_l(F)|');
    axis tight;
    grid on;
end

% find f' where |Z(F)| is maximized
[max_Z pos_max_Z] = max(Z);
f_star = f_axis(pos_max_Z);

% f_star = -m*Df = -m*DF*T => DF = -f_star/(m*T)
calculated_freq_offset = -f_star/(train_step*T);

%% cancel out the frequency offset
corr_exp = exp(-j*2*pi*f_star*[0:N-1]/train_step);   % "correction" exponent
r_bar = r_k.*corr_exp;    % multiply by the "correction" exponent

%verify cancelling of freq offset
if display_figures
    % plot real part of corrected symbols r_bar
    figure;
    subplot(2, 1, 1);
    stem(real(r_bar), 'r');
    title('In-phase CFO-corrected symbols');
    % plot imaginary part of corrected symbols r_bar
    subplot(2, 1, 2);
    stem(imag(r_bar), 'r');
    title('Quadrature CFO-corrected symbols');
    
    % scatterplot corrected symbols
    scatterplot(r_bar);
end

%% channel estimation and phase correction
tr = train_seq.';
h_hat = ((tr'*tr)^(-1))*(tr'*r_bar(tr_ind).');

r_bar_corr = conj(h_hat)*r_bar/abs(h_hat)^2;

%verify correction of phase
if display_figures
    % plot real part of corrected symbols r_bar_corr
    figure;
    subplot(2, 1, 1);
    stem(real(r_bar_corr), 'r');
    title('In-phase corrected symbols');
    % plot imaginary part of corrected symbols r_bar_corr
    subplot(2, 1, 2);
    stem(imag(r_bar_corr), 'r');
    title('Quadrature corrected symbols');
    
    % scatterplot corrected symbols
    scatterplot(r_bar_corr);     % now we can detect successfully!! :)
end

%% detection and errors measurement
s_k_I = sign(real(r_bar_corr));
s_k_Q = sign(imag(r_bar_corr));

s_k = s_k_I + j*s_k_Q;

errors = sum(a_n ~= s_k)



