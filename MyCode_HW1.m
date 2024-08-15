clc; clear; close all

% Define parameters
parameter.N_value = linspace(20, 200, 10); % Number of variables
parameter.rhos = linspace(0, 0.99, 5); % Correlation coefficients

% parameter.N_value = linspace(20, 200, 10); % Number of variables
% parameter.rhos = 0.7; % Correlation coefficients


parameter.MSE = zeros(length(parameter.N_value), length(parameter.rhos)); % MSE values

MonteC = 500; % Number of Monte Carlo iteration
% Loop over the number of variables and Monte Carlo simulations
for n = 1:length(parameter.N_value)
    N = parameter.N_value(n);
    for r = 1:length(parameter.rhos)
        rho = parameter.rhos(r);
        sample_cov_matrices = zeros(N, N, MonteC);
        Sigma = rho .^ abs((1:N)' - (1:N)); % True covariance matrix
        
        for i = 1:MonteC
            avg = zeros(N, 1);  
            X = mvnrnd(avg, Sigma, 1); % 1 means we generate only one sample at a time in each Monte-Carlo iteration. 
                                       % Therefore, we are in a single
                                       % sample oburvation. 

            sample_cov_matrices(:, :, i) = (X' * X) ./ (N-1); % sample covatiance matrix calculation
        end
        
        avg_sample_cov_matrix = mean(sample_cov_matrices, 3); % Average sample covariance matrix
        parameter.MSE(n, r) = norm(avg_sample_cov_matrix - Sigma, 'fro' )^2 / numel(Sigma);
    end
end

%%
% Plot MSE vs rho
figure;
for n = 1:length(parameter.N_value)
    semilogy(parameter.rhos, parameter.MSE(n, :), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'auto');
    hold on;
end

xlabel('\rho', 'FontWeight', 'bold');
ylabel('MSE', 'FontWeight', 'bold');
legend(arrayfun(@(x) ['N = ', num2str(x)], parameter.N_value, 'UniformOutput', false), 'Location', 'best', 'FontWeight', 'bold');
title('Mean Squared Error vs Correlation Coefficient', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

%%

% Plot MSE vs N
figure;
for r = 1:length(parameter.rhos)
    semilogy(parameter.N_value, parameter.MSE(:, r), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'auto');
    hold on;
end
xlabel('N', 'FontWeight', 'bold');
ylabel('MSE', 'FontWeight', 'bold');
legend(arrayfun(@(x) ['\rho = ', num2str(x)], parameter.rhos, 'UniformOutput', false), 'Location', 'best', 'FontWeight', 'bold');
title('Mean Squared Error vs Number of Measurements', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

%% ====================== 1.2.1: Frequency Estimation (L=1) =================
%  % Frequency Estimation for L=1
clear variables
close all
clc
Approach = 1;       %   1: Cww1   ||   2: Cww2: rho=0.9
% Parameters
L = 1;  % Number of sinusoids
rho=0.9;
N = [128 512 1024];  % Number of samples in each measurement
NMC = 500;  % Number of Monte Carlo simulations
n_us = 8;  % Up-sampling factor
a = 1;  % Amplitude
SNR = -20:2:30;  % SNR in dB scale
NSNR = length(SNR);  % Number of SNR values

% Function definitions
waveform = @(a,w,N_vec,phi) a*cos(w*(0:N_vec-1) + phi);  % Sine waveform
RandL = @() 2*pi*rand;  % Random phase generation

% Pre-allocate arrays for MSE and CRB
all_MSE = zeros(length(N), length(pi*[5/N(1), 10/N(1), 1/4, 1/2]), NSNR);
all_CRB = zeros(length(N), length(pi*[5/N(1), 10/N(1), 1/4, 1/2]), NSNR);

% Loop over different N values
for n = 1:length(N)
    N_i = N(n);
    w = pi*[5/N_i, 10/N_i, 1/sqrt(5), 1/2];  % Frequencies to estimate
    for i_w = 1:length(w)
        w_try = w(i_w);

        % Frequency axis for power spectral density
        f_ax = (-N_i/2:1/n_us:N_i/2-1/n_us)*2*pi/N_i;  % Upsampled frequency axis

        % Try different SNR values
        for nS = 1:NSNR
            sigma_w = (a^2) / (10^(SNR(nS)/10));  % Noise variance


            if Approach == 1
            % if Approach 1: Cww1
            C_ww = sigma_w .* eye(N_i);  % Covariance matrix (identity)
            else
            % if Approach 2: Cww2
            C_ww = sigma_w*(rho .^ abs((1:N)' - (1:N))); % True covariance matrix
            end


            MSE = zeros(NMC, 1);  % Initialize MSE matrix for current SNR

            % Monte Carlo simulations
            for nmc = 1:NMC
                phi_try = RandL();  % Generate random phase
                x = waveform(a, w_try, N_i, phi_try) + sqrt(sigma_w) * randn(1, N_i); % generating the waveform and sum the noise

                % Power Spectrum Density (PSD) Computation
                PSDx = fftshift((2/N_i) .* abs(fft(x, n_us*N_i)).^2);

                % Positive part of the PSD
                posPSD = PSDx((N_i*n_us)/2:end);
                pos_f_ax = f_ax((N_i*n_us)/2:end);
                [S0, maxIdx] = max(posPSD);

                w0_est = pos_f_ax(maxIdx);  % Estimated frequency
                MSE(nmc) = (w_try - w0_est)^2 / w_try^2;  % Correct MSE calculation (normalized)
            end

            all_MSE(n, i_w, nS) = mean(MSE);  % Store mean MSE for current settings
            all_CRB(n, i_w, nS) = ((1/2) * (1 / 10^(SNR(nS)/10))) / (N_i * (N_i - 1) * N_i);  % Cramer-Rao Bound (CRB)
        end
    end
end



%% Plotting
close all

% Plot MSE vs SNR for all N and w values
figure;


vertical_ax=1; % means which parameter to plot on vertical axis
% 1 ==> MSE
% 2 ==> CRB


plot_type = 1;
% 1 ==> MSE/CRB vs SNR for all N 
% 2 ==> MSE/CRB vs SNR for all w


iiw=pi*[5/N(n), 10/N(n), 1/4, 1/2];



if plot_type == 1  % plot type == 2 : for all N
for n = 1:length(N)

if vertical_ax == 1
        semilogy(SNR, squeeze(all_MSE(n, 3, :)), 'LineWidth', 2); %  change second argument
else
        semilogy(SNR, squeeze(all_CRB(n, 4, :)), 'LineWidth', 2); %  change second argument
end



        hold on;
        % the second argument of all_MSE/CRB should be as below: 
        % w=1 ==> 5/N(n)
        % w=2 ==> 10/N(n)
        % w=3 ==> 1/sqrt(5)
        % w=4 ==> 1/2

xlabel('SNR (dB)', 'FontWeight', 'bold');

if vertical_ax == 1
ylabel('MSE', 'FontWeight', 'bold');
%title('Mean Squared Error vs SNR for all N and \omega=π/2', 'FontWeight', 'bold');
else
 ylabel('CRB', 'FontWeight', 'bold');
%title('Cramér-Rao Bound vs SNR for all N and \omega=π/2', 'FontWeight', 'bold');
end


legend(arrayfun(@(x) ['N = ', num2str(x)], N, 'UniformOutput', false), 'Location', 'best', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;
end

else % plot type == 2 : for all w

for ii_w = 1:length(w)

if vertical_ax == 1
        semilogy(SNR, squeeze(all_MSE(3, ii_w, :)), 'LineWidth', 2); % try to change first argument
else
   
        semilogy(SNR, squeeze(all_CRB(3, ii_w, :)), 'LineWidth', 2); % try to change first argument
      
end


        hold on;
        % the first argument of all_MSE/CRB should be as below: 
        % n = 1 ==> N = 120
        % n = 2 ==> N = 512
        % n = 3 ==> N = 1024

xlabel('SNR (dB)', 'FontWeight', 'bold');


if vertical_ax == 1
ylabel('MSE', 'FontWeight', 'bold');
title('Mean Squared Error vs SNR for all \omega and N=1024', 'FontWeight', 'bold');
else
 ylabel('CRB', 'FontWeight', 'bold');
%title('Cramér-Rao Bound vs SNR for all \omega and N=100', 'FontWeight', 'bold');
title('Cramér-Rao Bound vs SNR for  N=100', 'FontWeight', 'bold');
end

%legend(arrayfun(@(x) [' = ', num2str(x)], w, 'UniformOutput', false), 'Location', 'best', 'FontWeight', 'bold');

legend("\omega = 5π/N","\omega = 10π/N","\omega = π/4","\omega = π/2", 'Location', 'best', 'FontWeight', 'bold');

%legend("N = 100","N = 500","N = 1024", 'Location', 'best', 'FontWeight', 'bold');

set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

end
end



%% 
% Plot PSD vs frequency 
figure;
plot(f_ax, PSDx, 'LineWidth', 2);
%plot(pos_f_ax, posPSD, 'LineWidth', 2);
xlabel('Frequency (rad/sample)', 'FontWeight', 'bold');
ylabel('Power Spectral Density', 'FontWeight', 'bold');
title('Power Spectral Density vs Frequency', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;
%xlim([-pi pi]);
xlim([1.54 1.6]);

%% ====================== 1.2.2: Frequency Estimation (L=2) =================


