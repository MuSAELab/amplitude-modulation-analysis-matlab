function [x, x_epoched] = istrfft_spectrogram(spectrogram_data)
%[x, x_epoched] = istrfft_spectrogram(spectrogram_data)
%
%     Compute the inverse STFT spectrogram for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     spectrogram_data : Structure with STFT spectrogram data, created with strfft_spectrogram()
% 
%     Returns
%     -------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%     x_epoched   = Segments form the signal or set of signals utilized to
%                   create the spectrogram in spectrogram_struct
%
% Example:
% xi = randn(100,1);
% xi_stft = strfft_spectrogram(xi, 10, 10, 5);
% xo = istrfft_spectrogram(xi_stft);

% Load data from Spectrogram structure
rFFT_data = spectrogram_data.rFFT_spectrogram;
win_size = spectrogram_data.win_size_samples;
win_shift = spectrogram_data.win_shift_samples;

% Generate psd_struct, to use irfft_psd()
psd_struct.fs   = spectrogram_data.fs;
psd_struct.channel_names = spectrogram_data.channel_names;
psd_struct.freq_axis = spectrogram_data.freq_axis;
psd_struct.win_function = spectrogram_data.win_function;
psd_struct.n_samples = win_size;

% Initialize rFFT_slice and x_epoched variables
[n_windows, n_freqs, n_channels] =  size(rFFT_data);
rfft_slide = zeros(n_freqs ,n_channels);
x_epoched = zeros(win_size, n_channels ,n_windows);

for i_window = 1 : n_windows
    % rFFT slice from spectrogram
    rfft_slide(:,:) = rFFT_data(i_window, :, :);
    % Generate psd_struct, to use irfft_psd()
    psd_struct.rFFT = rfft_slide; 
    % ifft_psd from the rFFT data recovers the signal or set of signals 'x'
    x_tmp = irfft_psd(psd_struct);
    x_epoched(:, :, i_window) = x_tmp;
end
% Merge epoched data
x = iepoching(x_epoched, win_shift);




