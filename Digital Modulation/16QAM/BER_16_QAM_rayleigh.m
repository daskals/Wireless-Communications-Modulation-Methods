% Clear all the previously used variables and close all figures
clear all;
close all;
 
format long;

% Frame Length 'Should be multiple of four or else padding is needed'
bit_count = 1*1000;

% Range of SNR over which to simulate 
Eb_No = 0: 1: 10;

% Convert Eb/No values to channel SNR

SNR = Eb_No + 10*log10(4);

% Start the main calculation loop
for aa = 1: 1: length(SNR)
    
    % Initiate variables
    T_Errors = 0;
    T_bits = 0;
    
    % Keep going until you get 100 errors
    while T_Errors < 100
    
        % Generate some random bits
        uncoded_bits  = round(rand(1,bit_count));

        % Split the stream into 4 substreams
        B = reshape(uncoded_bits,4,length(uncoded_bits)/4);
        B1 = B(1,:);
        B2 = B(2,:);
        B3 = B(3,:);
        B4 = B(4,:);
        
        % 16-QAM modulator
        % normalizing factor
        a = sqrt(1/10);

        % bit mapping
        tx = a*(-2*(B3-0.5).*(3-2*B4)-j*2*(B1-0.5).*(3-2*B2));
        
         % Variance = 0.5 - Tracks theoritical PDF closely
        ray = sqrt((1/2)*((randn(1,length(tx))).^2+(randn(1,length(tx))).^2));
        
        % Include The Fading
        rx = tx.*ray;
        
        % Noise variance
        N0 = 1/10^(SNR(aa)/10);

        % Send over Gaussian Link to the receiver
        rx = rx + sqrt(N0/2)*(randn(1,length(tx))+1i*randn(1,length(tx)));
        
%---------------------------------------------------------------

        % Equaliser
        rx = rx./ray;
        
        % 16-QAM demodulator at the Receiver
        a = 1/sqrt(10);

        B5 = imag(rx)<0;
        B6 = (imag(rx)<2*a) & (imag(rx)>-2*a);
        B7 = real(rx)<0;
        B8 = (real(rx)<2*a) & (real(rx)>-2*a);
        
        % Merge into single stream again
        temp = [B5;B6;B7;B8];
        B_hat = reshape(temp,1,4*length(temp));
    
        % Calculate Bit Errors
        diff =  uncoded_bits - B_hat ;
        T_Errors = T_Errors + sum(abs(diff));
        T_bits = T_bits + length(uncoded_bits);
        
    end
    % Calculate Bit Error Rate
    BER(aa) = T_Errors / T_bits;
    disp(sprintf('bit error probability = %f',BER(aa)));
    
    % Plot the received Symbol Constellation
    figure;
    grid on;
    plot(rx,'x');
    xlabel('Inphase Component');
    ylabel('Quadrature Component');
    Title(['Constellation of Transmitted Symbols for SNR =',num2str(SNR(aa))]);


end
  
%------------------------------------------------------------

% Finally plot the BER Vs. SNR(dB) Curve on logarithmic scale 
% BER through Simulation

figure(1);
semilogy(Eb_No,BER,'xr-','Linewidth',2);
hold on;
xlabel('E_b / N_o (dB)');
ylabel('BER');
title('E_b / N_o Vs BER plot for 16-QAM Modualtion in Rayleigh Channel');

% Theoretical BER
figure(1);
theoryBerAWGN = 0.5.*erfc(sqrt((10.^(Eb_No/10))));
semilogy(Eb_No,theoryBerAWGN,'g-+','Linewidth',2);
grid on;
legend('Rayleigh', 'AWGN');
axis([Eb_No(1,1) Eb_No(end) 0.00001 1]);