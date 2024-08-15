% Parameters
a = 1; % Amplitude
phi = 0; % Initial phase
N = 22000; % Number of samples
fs = 22000; % Sampling frequency
n = 0:N-1; % Time vector
SNR_dB_1 = 5; % Low SNR in dB
SNR_dB_2 = 25; % High SNR in dB

% Calculate gamma based on the given conditions
gamma = 3 * pi / (8 * N);

% Calculate Noise Variance for low SNR
sigma_w_squared_1 = (a^2) / (10^(SNR_dB_1 / 10));
sigma_w_1 = sqrt(sigma_w_squared_1);

% Generate the noise for low SNR
w_1 = sigma_w_1 * randn(size(n));

% Generate the frequency modulated sinusoid for low SNR
x_1 = a * cos(gamma * n.^2 + phi) + w_1;

% Calculate Noise Variance for high SNR
sigma_w_squared_2 = (a^2) / (10^(SNR_dB_2 / 10));
sigma_w_2 = sqrt(sigma_w_squared_2);

% Generate the noise for high SNR
w_2 = sigma_w_2 * randn(size(n));

% Generate the frequency modulated sinusoid for high SNR
x_2 = a * cos(gamma * n.^2 + phi) + w_2;

% Play the sound for low SNR
sound(x_1, fs);
pause(length(x_1) / fs + 2); % Pause to allow sound to finish playing

% Play the sound for high SNR
sound(x_2, fs);
pause(length(x_2) / fs + 2 ); % Pause to allow sound to finish playing

% Plot the signals
figure;
subplot(2,1,1);
plot(n, x_1);
title(['Frequency Modulated Sinusoid with SNR = ', num2str(SNR_dB_1), ' dB']);
xlabel('Sample number n');
ylabel('x[n]');

subplot(2,1,2);
plot(n, x_2);
title(['Frequency Modulated Sinusoid with SNR = ', num2str(SNR_dB_2), ' dB']);
xlabel('Sample number n');
ylabel('x[n]');

% Plot the instantaneous frequency
instantaneous_frequency = 2 * gamma * n;

figure;
subplot(2,1,1);
plot(n, instantaneous_frequency);
title(['Instantaneous Frequency with SNR = ', num2str(SNR_dB_1), ' dB']);
xlabel('Sample number n');
ylabel('Frequency (Hz)');

subplot(2,1,2);
plot(n, instantaneous_frequency);
title(['Instantaneous Frequency with SNR = ', num2str(SNR_dB_2), ' dB']);
xlabel('Sample number n');
ylabel('Frequency (Hz)');
%% PSD Plot:

% Compute and plot the PSD of x_1 and x_2
figure;
subplot(2,1,1);
pwelch(x_1, [], [], [], fs);
title(['Power Spectral Density with SNR = ', num2str(SNR_dB_1), ' dB']);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

subplot(2,1,2);
pwelch(x_2, [], [], [], fs);
title(['Power Spectral Density with SNR = ', num2str(SNR_dB_2), ' dB']);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%% Play Sounds

% Play the sound for low SNR
sound(x_1, fs);
pause(length(x_1) / fs + 2); % Pause to allow sound to finish playing

% Play the sound for high SNR
sound(x_2, fs);
pause(length(x_2) / fs + 2 ); % Pause to allow sound to finish playing