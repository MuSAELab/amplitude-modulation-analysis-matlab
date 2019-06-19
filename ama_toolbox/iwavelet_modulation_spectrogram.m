function x = iwavelet_modulation_spectrogram(modulation_spectrogram_data)
% x = istrfft_modspectrogram(modulation_spectrogram_data)
%
%     Compute the inverse CWT-based modulation spectrogram for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     modulation_spectrogram_data : Structure with CWT-based modulation spectrogram data, 
%           created with wavelet_modulation_spectrogram()
% 
%     Returns
%     -------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%
% Example:
% xi = randn(1000,1);
% xi_cwt = wavelet_modulation_spectrogram(xi, 1000);
% xo = iwavelet_modulation_spectrogram(xi_cwt);

% Number of channels from Modspectrogram structure
n_channels = size(modulation_spectrogram_data.rFFT_modulation_spectrogram, 3);

% Prepare psd_tmp_struct to perform irFFT on Modulation Spectogram
psd_tmp_data.freq_axis = modulation_spectrogram_data.freq_mod_axis;
psd_tmp_data.fs = modulation_spectrogram_data.fs_mod;
psd_tmp_data.win_function = modulation_spectrogram_data.win_function_x;
psd_tmp_data.n_samples = modulation_spectrogram_data.n_samples;


for i_channel = 1 : n_channels
    % Slide with the rFFT coeffients of the 2nd FFT 
    psd_tmp_data.rFFT = transpose(modulation_spectrogram_data.rFFT_modulation_spectrogram(:,:,i_channel));   
    % Recovers the Square Root of the Power Spectrogram
    sqrt_pwr_spectrogram = irfft_psd(psd_tmp_data);
    
    % Recovers the Magnitude of the Wavelet Coefficients
    pwr_spectrogram = sqrt_pwr_spectrogram .^2;
    pwr_spectrogram = pwr_spectrogram * modulation_spectrogram_data.fs_mod *  modulation_spectrogram_data.n_samples;
    pwr_spectrogram = pwr_spectrogram ./2 ;
    spectrogram_abs = sqrt(pwr_spectrogram);
        
    % Recovers the Angle values of the Spectrogram
    spectrogram_angle = angle(modulation_spectrogram_data.spectrogram_data.wavelet_coefficients(:,:,i_channel));
    
    % Creates the rFFT coefficients of the 1st FFTs
    modulation_spectrogram_data.spectrogram_data.wavelet_coefficients(:,:,i_channel) = spectrogram_abs .* exp(1i .* spectrogram_angle );   
end

% Recovers the origial signal or set of signals
x = iwavelet_spectrogram(modulation_spectrogram_data.spectrogram_data);

