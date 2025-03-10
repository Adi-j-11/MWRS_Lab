nFFT = 64; % FFT size
nDSC = 52; % Number of data subcarriers
nSym = 1e4; % Number of symbols
EbN0dB = 0:10; % Bit-to-noise ratio in dB
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % Convert to symbol SNR

nErr = zeros(1, length(EbN0dB)); % Initialize error storage

for ii = 1:length(EbN0dB)
    % Generate random bit stream
    ipBit = randi([0,1], 1, nDSC * nSym);  
    
    % BPSK Mapping (0 → -1, 1 → +1)
    ipMod = 2 * ipBit - 1;  

    % Reshape into OFDM symbols
    ipMod = reshape(ipMod, nDSC, nSym).';  

    % Assign to subcarriers (-26 to -1, +1 to +26)
    xF = [zeros(nSym,6), ipMod(:,1:nDSC/2), zeros(nSym,1), ipMod(:,nDSC/2+1:nDSC), zeros(nSym,5)];  

    % IFFT and normalize power
    xt = (nFFT / sqrt(nDSC)) * ifft(fftshift(xF.')).';  

    % Add cyclic prefix (last 16 samples prepended)
    xt = [xt(:,49:64), xt];  

    % Convert to serial stream
    xt = reshape(xt.', 1, []);  

    % Generate complex Gaussian noise
    nt = (randn(1, length(xt)) + 1j * randn(1, length(xt))) / sqrt(2);  

    % Add noise (accounting for cyclic prefix power loss)
    yt = sqrt(80/64) * xt + 10^(-EsN0dB(ii)/20) * nt;  

    % Reshape received signal into symbols
    yt = reshape(yt, 80, []).';  

    % Remove cyclic prefix
    yt = yt(:, 17:80);  

    % FFT and shift
    yF = (sqrt(nDSC) / nFFT) * fftshift(fft(yt.')).';  

    % Extract data subcarriers
    yMod = yF(:, [7:32, 34:59]);  

    % BPSK Demodulation (Threshold at 0)
    ipBitHat = real(yMod) > 0;  

    % Reshape back to bit sequence
    ipBitHat = reshape(ipBitHat.', 1, []);  

    % Count bit errors
    nErr(ii) = sum(ipBitHat ~= ipBit);  
end

% Compute simulated BER
simBer = nErr / (nDSC * nSym);  

% Compute theoretical BER for BPSK
theoryBer = (1/2) * erfc(sqrt(10.^(EbN0dB / 10)));  

% Plot results
figure;  
semilogy(EbN0dB, theoryBer, 'bs-', 'LineWidth', 2);  
hold on;  
semilogy(EbN0dB, simBer, 'mx-', 'LineWidth', 2);  
grid on;  

xlabel('E_b/N_0 (dB)');  
ylabel('Bit Error Rate');  
legend('Theory', 'Simulation');  
title('BER for BPSK using OFDM');  
