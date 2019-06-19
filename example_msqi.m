%% Example MSQI
% This example shows the use of the amplitude modulation analysis toolkit
% to compute the Modulation Spectrum-Based ECG Quality Index (MSQI) as a 
% blind metric to measure the signal-to-noise (SNR) ration in ECG signals
% 
% The MSQI was originally presented in:
%
% D. P. Tobon V., T. H. Falk, and M. Maier, "MS-QI:  A  Modulation
% Spectrum-Based ECG Quality Index for Telehealth Applications", IEEE
% Transactions on Biomedical Engineering, vol. 63, no. 8, pp. 1613-1622,
% Aug. 2016

%% ECG signal 
% load ECG signal
% x_clean is a 5-s segment of clean ECG signal
% x_noisy is the x_clean signal contaminated with pink noise to have a 0db SNR
load('./example_data/msqi_ecg_data.mat');

%% Compute MSQI and heart rate (HR), and plot modulation spectrogram
[msqi_clean, hr_clean, modulation_spectrogram_clean] = msqi_ama(x_clean, fs);
[msqi_noisy, hr_noisy, modulation_spectrogram_noisy] = msqi_ama(x_noisy, fs);

fprintf('HR = %f bpm \r\n', hr_clean);
fprintf('MSQI for clean ECG = %f \r\n', msqi_clean);
fprintf('MSQI for noisy ECG = %f \r\n', msqi_noisy);

%% Plot modulation spectrograms      
figure('units','normalized','outerposition',[0 0 1 1])
% clean 
subplot(3,2,1)
plot_signal(x_clean, fs)
title('Clean ECG');
subplot(3,2,[3,5])
plot_modulation_spectrogram_data(modulation_spectrogram_clean, ...
    [], [0, 60], [], [-90, -40])
title('Modulation spectrogram clean ECG');
axis square
subplot(3,2,2)
plot_signal(x_noisy, fs)
title('Noisy ECG');
subplot(3,2,[4,6])
plot_modulation_spectrogram_data(modulation_spectrogram_noisy, ...
    [], [0, 60], [], [-90, -40])
title('Modulation spectrogram noisy ECG');
axis square
