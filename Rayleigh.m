clc; clear; close all;

% Simulation parameters
N = 20000;  % Number of samples
SNR_limit = 35;
SNR_db = -5:0.5:SNR_limit; 
SNR = 10.^(SNR_db/10); 
var = 1;  % Rayleigh variance (sigma^2)
nstd = sqrt(var);

% Generate Rayleigh random variable
u = rand(1, N);
r = sqrt(-2 * var * log(u)); 

% Plot histogram of Rayleigh random variable
figure;
histogram(r, 100, 'Normalization', 'pdf'); 
title('Rayleigh Random Variable Histogram');
xlabel('Random Variable R');
ylabel('Frequency');

% Rayleigh Probability Density Function (PDF)
a = 0:0.01:10;
R = (a / var) .* exp(-(a.^2) / (2 * var));

figure;
plot(a, R, 'r', 'LineWidth', 2);
title('Rayleigh PDF');
xlabel('Random Variable');
ylabel('Probability');
legend('Variance = 1');

% Compute Theoretical BER for fading and AWGN
Pe_BPSK_fading = 0.5 * (1 - sqrt((var * SNR) ./ (1 + var * SNR)));  
Pe_BFSK_fading = 0.5 * (1 - sqrt(var * SNR ./ (2 + var * SNR)));  
Pe_DPSK_fading = 0.5 ./ (1 + var * SNR);

Pe_BPSK_AWGN = 0.5 * erfc(sqrt(SNR));  
Pe_BFSK_AWGN = 0.5 * erfc(sqrt(SNR / 2));  
Pe_DPSK_AWGN = 0.5 * exp(-SNR);  

% Plot BER performance
figure;
semilogy(SNR_db, Pe_BPSK_fading, 'r.-', ...
         SNR_db, Pe_BFSK_fading, 'r*-', ...
         SNR_db, Pe_DPSK_fading, 'r--', ...
         SNR_db, Pe_BPSK_AWGN, 'b.-', ...
         SNR_db, Pe_BFSK_AWGN, 'b*-', ...
         SNR_db, Pe_DPSK_AWGN, 'b--', 'LineWidth', 1.5);

axis([-5 SNR_limit 1e-6 1]);
grid on;
title('Performance of BPSK, BFSK, and DPSK');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('BPSK (Fading)', 'BFSK (Fading)', 'DPSK (Fading)', ...
       'BPSK (AWGN)', 'BFSK (AWGN)', 'DPSK (AWGN)');

% Monte Carlo Simulation for BPSK in Rayleigh fading
Eb = 1; 
EbNo_dB = 0:5:35; 
No_over_2 = Eb * 10.^(-EbNo_dB/10);
BER = zeros(1, length(EbNo_dB));

for i = 1:length(EbNo_dB)
    no_errors = 0;
    no_bits = 0;

    while no_errors <= 10  
        u = rand;
        alpha = sqrt(-2 * var * log(u));  
        noise = sqrt(No_over_2(i)) * randn;  
        y = alpha * sqrt(Eb) + noise;  
        y_d = y <= 0;  
        
        no_bits = no_bits + 1;
        no_errors = no_errors + y_d;
    end
    BER(i) = no_errors / no_bits;  
end

% Theoretical BER for BPSK in Rayleigh fading
rho_b = Eb ./ No_over_2 * var;
P2 = 0.5 * (1 - sqrt(rho_b ./ (1 + rho_b)));

% Plot Monte Carlo vs Theoretical BER
figure;
semilogy(EbNo_dB, BER, '-*', EbNo_dB, P2, '-o', 'LineWidth', 1.5);
grid on;
title('Monte Carlo Simulation for BPSK in Rayleigh Fading');
xlabel('SNR per Bit (dB)');
ylabel('Bit Error Rate (BER)');
legend('Monte Carlo Simulation', 'Theoretical Value');
