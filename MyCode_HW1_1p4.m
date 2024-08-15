load('Hw2a.mat');
x=Sin_d;
%%
% Play the sound 
sound(x, fs);
pause(length(x) / fs + 2); % Pause to allow sound to finish playing

%%
close all
% Compute and plot the PSD of x_a and x_b
figure;
pwelch(x, [], [], [], fs);
%pwelch(x);
title('Power Spectral Density');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
%%
spectrogram(x(:,1))

