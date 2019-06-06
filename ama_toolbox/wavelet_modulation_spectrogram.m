function modulation_spectrogram_data = wavelet_modulation_spectrogram( x, fs, n_cycles, freq_vct, fft_factor_x, win_funct_x, channel_names)
%modulation_spectrogram_data = WAVELET_MODSPECTROGRAM( x, fs, n_cycles, freq_vct, fft_factor_x, win_funct_x, channel_names)
%     Compute the Modulation Spectrogram using the Wavelet for one or a set of REAL signals 'x'.
%         
%     Parameters
%     ----------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%     fs : Sampling frequency 
%          in Hz
%     n : Number of cicles inside the Gaussian curve 
%         (Default 6)
%     freq_vct : 1D array 
%         with frequencies to compute the CWT (Default = [1 : 1 : fs/2] )
%     fft_factor_x : Number of elements to perform the FFT is given as:
%         n_fft_x  = fft_factor_x * n_samples, (default, fft_factor_x = 1)
%     win_function_x : Window to apply in the rFFT 
%         (Default 'Hamming')   
%     channel_names : Names of the signals
%         (Default Signal-XX with XX 1, 2, ... n_channels) 
% 
%     Returns
%     -------
%     modulation_spectrogram_data : Dictionary with Modulation Spectrogram data, with the elements:
%        rFFT_modulation_spectrogram
%            rFFT values for each window (u), scaled by the Window RMS       
%        power_modulation_spectrogram :
%            Power modulation spectrogram (u^2 / Hz) 
%        fs : 
%            Sampling frequency (Hz)
%        fs_mod : 
%            Sampling frequency of modulation-frequency (Hz)         
%        freq_axis :
%            Frequency axis for rFFT and PSD (Hz)
%        freq_delta :
%            Frequency axis step (Hz)
%        freq_mod_axis :
%            Modulation-frequency axis for rFFT_modspec and pwr_modspec (Hz)     
%        freq_mod_delta :
%            Modulation-frequency step (Hz)
%        n_fft_x :
%            Number of elements utilized to perform the FFT
%        win_funct_x :
%            Window to apply in the 2nd rFFT                      
%        n_samples :
%            Number of samples of the signal or signals 'x'
%        spectrogram_data : 
%            Dictionary with Spectrogram data
%        channel_names :
%            Names of channels

% get class of x
x_class = class(x);

% validate 'freq_vct' argument
if ~exist('freq_vct','var') || isempty(freq_vct)
    freq_vct = 1 : floor(fs / 2);
end

% validate 'n_cycles' argument
if ~exist('n_cycles','var') || isempty(n_cycles)
    n_cycles = 6;
end

% validate 'win_funct_x' argument
if ~exist('win_funct_x','var') || isempty(win_funct_x)
    win_funct_x = 'hamming';
end

% validate 'fft_factor_x' argument
if ~exist('fft_factor_x','var') || isempty(fft_factor_x)
    fft_factor_x = 1;
end

% validate 'channel_names' argument
if ~exist('channel_names','var') || isempty(channel_names)
    channel_names = {};
end

% compute wavelet spectrogram
spectrogram_data = wavelet_spectrogram(x, fs, n_cycles, freq_vct, channel_names);
[n_windows, n_freqs, n_channels] =  size(spectrogram_data.wavelet_coefficients);

% number of elements for FFT of the spectrogram
n_fft_x =  fft_factor_x * n_windows;

% generate default channel names
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('Signal-%02d',ic);
    end
end

% modulation sampling frequency
fs_mod = fs;

% the AM analysis is made in the Amplitude derived from the Power Spectrogram
for i_channel = 1 : n_channels
  % data to generate the Modulation Spectrogram
    spectrogram_1ch = sqrt(spectrogram_data.power_spectrogram(:,:,i_channel)) ;
    %spectrogram_1ch = abs(spectrogram_data.rFFT_spectrogram(:,:,i_channel));
    %spectrogram_1ch = angle(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));
    %spectrogram_1ch = imag(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));
    %spectrogram_1ch = real(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));

    % Compute rfft_psd on each frequency timeseries
    mod_psd_struct = rfft_psd(spectrogram_1ch, fs_mod, n_fft_x, win_funct_x, channel_names );

    if i_channel == 1
        % modulation frequency axis
        fmod_ax = mod_psd_struct.freq_axis;
        % modulation frequency delta
        fmod_delta = mod_psd_struct.freq_delta;

        % initialize 'rFFT_modspec'  and 'pwr_modspec'
        n_freqsmod = numel(fmod_ax);
        rFFT_modspec = zeros(n_freqs, n_freqsmod ,n_channels, x_class);
        pwr_modspec  = zeros(n_freqs, n_freqsmod ,n_channels, x_class);
    end
    % rFFT data
    rFFT_modspec(:, :, i_channel) = transpose(mod_psd_struct.rFFT);
    % power data
    pwr_modspec(:, :, i_channel) = transpose(mod_psd_struct.PSD);
end

% scale 'pwr_modspec' by modulation-frequency delta
pwr_modspec = pwr_modspec / fmod_delta;

% output 'modulation_spectrogram_data' structure
modulation_spectrogram_data.rFFT_modulation_spectrogram = rFFT_modspec;
modulation_spectrogram_data.power_modulation_spectrogram = pwr_modspec;
modulation_spectrogram_data.fs = fs;
modulation_spectrogram_data.fs_mod = fs_mod;
modulation_spectrogram_data.freq_axis = spectrogram_data.freq_axis;
modulation_spectrogram_data.freq_delta = spectrogram_data.freq_delta;
modulation_spectrogram_data.freq_mod_axis = fmod_ax;
modulation_spectrogram_data.freq_mod_delta = fmod_delta;
modulation_spectrogram_data.n_fft_x = n_fft_x;
modulation_spectrogram_data.win_funct_x = win_funct_x;
modulation_spectrogram_data.n_samples = spectrogram_data.n_samples;
modulation_spectrogram_data.spectrogram_data = spectrogram_data;
modulation_spectrogram_data.channel_names = channel_names;
end
