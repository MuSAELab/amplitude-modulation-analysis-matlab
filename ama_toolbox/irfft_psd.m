function x = irfft_psd(psd_data)
%x = irfft_psd(psd_data)
%
%     Compute the inverse PSD for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     psd_data : Structure with PSD data, created with rfft_psd()
% 
%     Returns
%     -------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%
% Example:
% xi = randn(100,1);
% xi_psd = rfft_psd(xi, 10);
% xo = irfft_psd(xi_psd);

% Load data from PSD structure
rFFT_data = psd_data.rFFT;
f_ax = psd_data.freq_axis;
fs = psd_data.fs;
win_function = psd_data.win_function;
n_samples = psd_data.n_samples;
n_channels = size(rFFT_data, 2);

% Find the number of elements used for the rFFT
if f_ax(end) < fs/2
    % elements for FFT was odd
    n_fft = (numel(f_ax) * 2) - 1;
elseif f_ax(end) - fs/2 < 1000 * eps
    % elements for FFT was even
    n_fft = (numel(f_ax) - 1) * 2;
end

% Window RMS
win = window(win_function, n_samples);
win_rms = sqrt(sum(win.^2) / n_samples);

% IRFFT
X = rFFT_data * win_rms;
x_tmp = irfft(X, n_fft);

% Keep only n_samples points
x = x_tmp(1 : n_samples, :);

% Un-Windowing
win = window(win_function,n_samples);
win_mat = repmat(win, 1, n_channels);
x = x ./ win_mat;
