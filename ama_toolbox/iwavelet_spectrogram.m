function x = iwavelet_spectrogram(spectrogram_data)
%[x, x_epoched] = iwavelet_spectrogram(spectrogram_data)
%
%     Compute the inverse CWT Spectrogram for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     spectrogram_data : Structure with CWT Spectrogram data, created with wavelet_spectrogram()
% 
%     Returns
%     -------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%
% Example:
% xi = randn(1000,1);
% xi_cwt = wavelet_spectrogram(xi, 1000);
% xo = iwavelet_spectrogram(xi_cwt);

% compute the scaling factor for each wavelet kernel
s = spectrogram_data.n_cycles ./ ( 2 * pi * spectrogram_data.freq_axis);
A = 1 ./ ((s.^2) * pi).^(1/4);

% compute the mean across scaled "filtered" signals
x = squeeze(mean( bsxfun(@rdivide, real(spectrogram_data.wavelet_coefficients) , A ), 2));




