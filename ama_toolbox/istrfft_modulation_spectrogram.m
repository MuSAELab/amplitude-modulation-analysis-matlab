function x = istrfft_modulation_spectrogram(modulation_spectrogram_data)
% x = istrfft_modulation_spectrogram(modulation_spectrogram_data)
%
%     Compute the inverse STFT-based modulation spectrogram for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     modulation_spectrogram_data : Structure with STFT-based modulation spectrogram data, 
%           created with strfft_modulation_spectrogram()
% 
%     Returns
%     -------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%
% Example:
% xi = randn(100,1);
% xi_mod_stft = strfft_modulation_spectrogram(xi, 10, 10, 5);
% xo = istrfft_modulation_spectrogram(xi_mod_stft);

% Number of channels from Modspectrogram structure
n_channels = size(modulation_spectrogram_data.rFFT_modulation_spectrogram, 3);

% Prepare psd_tmp_data to perform irFFT on Modulation Spectogram
psd_tmp_data.freq_axis = modulation_spectrogram_data.freq_mod_axis;
psd_tmp_data.fs = modulation_spectrogram_data.fs_mod;
psd_tmp_data.win_function = modulation_spectrogram_data.win_function_x;
psd_tmp_data.n_samples = modulation_spectrogram_data.n_windows;


for i_channel = 1 : n_channels
    % Slide with the rFFT coeffients of the 2nd FFT 
    psd_tmp_data.rFFT = transpose(modulation_spectrogram_data.rFFT_modulation_spectrogram(:,:,i_channel));   
    % Recovers the Square Root of the Power Spectrogram
    sqrt_pwr_spectrogram = irfft_psd(psd_tmp_data);
    % Power Spectrogram
    pwr_spectrogram = sqrt_pwr_spectrogram .^ 2;
    % Scale Power Spectrogram by (n_windows * time_delta)
    pwr_spectrogram = pwr_spectrogram * modulation_spectrogram_data.spectrogram_data.n_windows * modulation_spectrogram_data.spectrogram_data.time_delta;
    % Scale Power Spectrogram by (freq_delta)
    pwr_spectrogram = pwr_spectrogram * modulation_spectrogram_data.spectrogram_data.freq_delta;
    % Scale Power Spectrogram by the number of samples used
    pwr_spectrogram = pwr_spectrogram / (1 / modulation_spectrogram_data.spectrogram_data.n_fft .^2);
    % Divde by 2 all the elements except DC and the Nyquist point (in even case)  
    pwr_spectrogram = pwr_spectrogram / 2;
    pwr_spectrogram(:, 1) = pwr_spectrogram(:, 1) * 2;
    if ~mod(modulation_spectrogram_data.spectrogram_data.n_fft, 2)
        % NFFT was even, then 
        pwr_spectrogram(:, end) = pwr_spectrogram(:, end) * 2;
    end
    spectrogram_abs = sqrt(pwr_spectrogram);
    % Recovers the Angle values of the Spectrogram
    spectrogram_angle = angle(modulation_spectrogram_data.spectrogram_data.rFFT_spectrogram(:,:,i_channel));
    % Creates the rFFT coefficients of the 1st FFTs
    modulation_spectrogram_data.spectrogram_data.rFFT_spectrogram(:,:,i_channel) = spectrogram_abs .* exp(1i .* spectrogram_angle );   
end

% Recovers the origial signal or set of signals
x = istrfft_spectrogram(modulation_spectrogram_data.spectrogram_data);

