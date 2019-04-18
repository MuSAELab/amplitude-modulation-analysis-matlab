%% Example 04
% This example shows the amplitude modulation analysis toolbox for
% speech data
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

%% Speech signal
% The speech signal p234_004.wav is one sample from the: 
% CSTR VCTK Corpus: English Multi-speaker Corpus for CSTR Voice Cloning Toolkit
% avialable in: https://datashare.is.ed.ac.uk/handle/10283/2651
[x,fs] = audioread('./example_data/p234_004.wav');
x_name = 'speech';
% 1s segment to analyze
x = x((fs*1.6)+1 : fs*3.6);
sound(x,fs)

%% STFT-based
% Parameters
win_size_sec = 0.04;   % window length for the STFFT (seconds)
win_shft_sec = 0.01;  % shift between consecutive windows (seconds)

figure()
subplot(4,5,[1:5])
plot_signal(x, fs, x_name ); colorbar;

stft_spectrogram = strfft_spectrogram(x, fs, round(win_size_sec*fs), round(win_shft_sec*fs), [], [], {x_name});
subplot(4,5,[6:10])
plot_spectrogram_data(stft_spectrogram)

stft_modulation_spectrogram = strfft_modulation_spectrogram(x, fs, round(win_size_sec*fs), round(win_shft_sec*fs), [], [], 1, [], {x_name});
subplot(4,5,[12:14, 17:19])
plot_modulation_spectrogram_data(stft_modulation_spectrogram,1 ,[], [0,20],[-90, -50])


%% Parameters for CWT for speech signal 
n_cycles = 6;          % number of cycles (for Complex Morlet)
up_lim = floor(log(fs/2) / log(2));
frequency_vector = 2.^[1:0.2:up_lim]; % vector of frequencies to compute the CWT

figure()
subplot(4,5,[1:5])
plot_signal(x, fs, x_name ); colorbar;

cwt_spectrogram = wavelet_spectrogram(x, fs, n_cycles, frequency_vector, {x_name});
subplot(4,5,[6:10])
plot_spectrogram_data(cwt_spectrogram)

cwt_modulation_spectrogram = wavelet_modulation_spectrogram(x, fs, n_cycles, frequency_vector, [], [], {x_name});
subplot(4,5,[12:14, 17:19])
plot_modulation_spectrogram_data(cwt_modulation_spectrogram, 1, [], [0,20],[-90, -50])
