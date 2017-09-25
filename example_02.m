%% Example 02
% Script to show the use of the functions:
%
% rfft_psd()                        Compute PSD using rFFT
% strfft_spectrogram()              Compute Spectrogram using STFFT
% strfft_modulation_spectrogram()   Compute Modulation Spectrogram using STFFT
% wavelet_spectrogram()             Compute Spectrogram using wavelet transformation
% wavelet_modulation_spectrogram()  Compute Modulation Spectrogram using wavelet transformation
%
% plot_signal()                      Plot a signal in time domain
% plot_psd_data()                    Plot PSD data obtained with rfft_psd()
% plot_spectrogram_data()            Plot Spectrogram data obtained with
%                                         strfft_spectrogram() or wavelet_spectrogram()
% plot_modulation_spectrogram_data() Plot Modulation Spectrogram data obtained with
%                                         strfft_modspectrogram() or wavelet_modspectrogram()
%
% Moreover, this script compares diverse ways to compute the power of a
% signal in the Time, Frequency, Time-Frequency and Frequency-Frequency domains
%
close all
clear

%% Test signal

fs = 500;
T = 1/fs;
t1_vct = (0 : 1/fs : 10 - 1/fs)';

x1 = 3 * sin (2 * pi * 10 * t1_vct);
x2 = 2 * sin (2 * pi * 24 * t1_vct);
x3 = 1 * randn([numel(t1_vct), 1]);

x = [x1; x2; x3];
n = numel(x);

%% Plot signal
plot_signal(x, fs, 'test-signal');

%% Power in Time Domain
% Energy of the signal
energy_x = T .* sum(x.^2);
duration = T .* n;

% Power of the signal
power_x = energy_x ./ duration;

% A simpler way is
power_x_2 = (1 / n) * sum(x.^2);

%% Power in Frequency domain
% Power using FFT
X = fft(x);
power_X = (1 / n.^2) * sum(X.*conj(X));

% Power using its PSD from rFFT
psd_rfft_r = rfft_psd(x, fs, [], 'rectwin');
f_step = psd_rfft_r.freq_axis(2);
power_psd_rfft_x_rw = f_step * sum(psd_rfft_r.PSD);
figure()
plot_psd_data(psd_rfft_r)

% Power using its PSD from rFFT
psd_rfft_b = rfft_psd(x, fs, [], 'blackmanharris');
f_step = psd_rfft_b.freq_axis(2);
power_psd_rfft_x_bh = f_step * sum(psd_rfft_b.PSD);
figure()
plot_psd_data(psd_rfft_b)

%% Power from STFFT Spectrogram (Hamming window)
w_size =  1 * fs;
w_shift = 0.5 * w_size;
win_funct = 'hamming';
rfft_spect_h = strfft_spectrogram(x, fs, w_size, w_shift, [], win_funct );
power_spect_h = sum(sum(rfft_spect_h.power_spectrogram)) * rfft_spect_h.freq_delta * rfft_spect_h.time_delta;
figure()
plot_spectrogram_data(rfft_spect_h);

%% Power from STFFT Spectrogram (Rectangular window)
w_size =  1 * fs;
w_shift = 0.5 * w_size;
win_funct = 'rectwin';
rfft_spect_r = strfft_spectrogram(x, fs, w_size, w_shift, [], win_funct );
power_spect_r = sum(sum(rfft_spect_r.power_spectrogram)) * rfft_spect_r.freq_delta * rfft_spect_r.time_delta;
figure()
plot_spectrogram_data(rfft_spect_r);

%% Power from Wavelet Spectrogram N = 6
wav_spect_6 = wavelet_spectrogram(x, fs, 6);
power_wav_6 = sum(sum(wav_spect_6.power_spectrogram)) * wav_spect_6.freq_delta * wav_spect_6.time_delta ;
figure()
plot_spectrogram_data(wav_spect_6);

%% Power from Wavelet Spectrogram N = 10
wav_spect_10 = wavelet_spectrogram(x, fs, 10);
power_wav_10 = sum(sum(wav_spect_10.power_spectrogram)) * wav_spect_10.freq_delta * wav_spect_10.time_delta ;
figure()
plot_spectrogram_data(wav_spect_10);

%% Power from Modulation Spectrogram STFFT
w_size =  1 * fs;
w_shift = 0.5 * w_size;
rfft_mod_b = strfft_modulation_spectrogram(x, fs, w_size, w_shift, [], win_funct, [], win_funct);
power_mod = sum(sum(rfft_mod_b.power_modulation_spectrogram) * rfft_mod_b.freq_delta * rfft_mod_b.freq_mod_delta);
figure()
plot_modulation_spectrogram_data(rfft_mod_b)

%% Power from Modulation Spectrogram Wavelet
wav_mod_6 = wavelet_modulation_spectrogram(x, fs, 6, [], [], win_funct);
power_mod_w = sum(sum(wav_mod_6.power_modulation_spectrogram) * wav_mod_6.freq_delta * wav_mod_6.freq_mod_delta);
figure()
plot_modulation_spectrogram_data(wav_mod_6)
