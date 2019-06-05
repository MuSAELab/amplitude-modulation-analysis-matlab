function modulation_spectrogram_data = strfft_modulation_spectrogram( x, fs, win_size, win_shift, fft_factor_y, win_funct_y, fft_factor_x, win_funct_x, channel_names)
%modspect_struct = STRFFT_MODSPECTROGRAM( x, fs, win_size, win_shift, fft_factor_y, win_funct_y, fft_factor_x, win_funct_x, channel_names)
% Compute the Modulation Spectrogram for one or a set of REAL signals 'x'.
%         
%     Parameters
%     ----------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%     fs : Sampling frequency 
%          in Hz
%     win_size :
%         Size of the sliding window for STFFF (samples)
%     win_shift :
%         Shift between consecutive windows (samples)   
%     fft_factor_y : Number of elements to perform the 1st FFT is given as:
%         n_fft_y  = fft_factor_y * n_samples, (default, fft_factor_y = 1)
%     win_function_y : Window to apply in the 1st FFT 
%         (Default 'Hamming')
%     fft_factor_x : Number of elements to perform the 2nd FFT is given as:
%         n_fft_x  = fft_factor_x * n_samples, (default, fft_factor_x = 1)
%     win_function_x : Window to apply in the 2nd rFFT 
%         (Default 'Hamming')   
%     n_fft : Number of samples to compute the FFT
%             (Default = n_samples in array x)   
%     win_funct : Window function applied to the signal 
%         (Default 'Hamming')
%     channel_names : Names of the signals
%         (Default Signal-XX with XX 1, 2, ... n_channels) 
% 
%     Returns
%     -------
%     modulation_spectrogram_data : Structure with Modulation Spectrogram data, with the elements:
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
%        win_size_samples :
%            Size of the sliding window for STFFF (samples)
%        win_shift_samples :
%            Shift between consecutive windows (samples)   
%        n_fft_y :
%            Number of elements utilized to perform the 1st FFT
%        n_fft_x :
%            Number of elements utilized to perform the 2nd FFT
%        win_funct_y :
%            Window to apply in the 1st rFFT            
%        win_funct_x :
%            Window to apply in the 2nd rFFT                      
%        n_windows :
%            Number of ST windows
%        n_samples :
%            Number of samples of the signal or signals 'x'
%        spectrogram_data : 
%            Structure with Spectrogram data
%        channel_names :
%            Names of channels

% get class of x
x_class = class(x);

% validate 'win_funct_y' argument
if ~exist('win_funct_y','var') || isempty(win_funct_y)
    win_funct_y = 'hamming';
end

% validate 'win_funct_x' argument
if ~exist('win_funct_x','var') || isempty(win_funct_x)
    win_funct_x = 'hamming';
end

% validate 'fft_factor_y' argument
if ~exist('fft_factor_y','var') || isempty(fft_factor_y)
    fft_factor_y = 1;
end

% number of elements for the 1st FFT
n_fft_y = fft_factor_y * win_size;

% validate 'fft_factor_x' argument
if ~exist('fft_factor_x','var') || isempty(fft_factor_x)
    fft_factor_x = 1;
end

% validate 'channel_names' argument
if ~exist('channel_names','var') || isempty(channel_names)
    channel_names = {};
end

% compute STFFT spectrogram
spectrogram_data = strfft_spectrogram(x, fs, win_size, win_shift, n_fft_y, win_funct_y, channel_names);
[n_windows, n_freqs, n_channels] =  size(spectrogram_data.rFFT_spectrogram);
% Number of elements for the 2nd FFT
n_fft_x =  fft_factor_x * n_windows;

% Generate default channel names, if needed
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('Signal-%02d',ic);
    end
end

% modulation sampling frequency
fs_mod = 1 / (win_shift / fs);

% the AM analysis is made in the Amplitude derived from the Power Spectrogram
for i_channel = 1 : n_channels
    % data to generate the Modulation Spectrogram
    spectrogram_1ch = sqrt(spectrogram_data.power_spectrogram(:,:,i_channel)) ;
    %spectrogram_1ch = abs(spectrogram_data.rFFT_spectrogram(:,:,i_channel));
    %spectrogram_1ch = angle(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));
    %spectrogram_1ch = imag(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));
    %spectrogram_1ch = real(squeeze(spectrogram_data.rFFT_spectrogram(:,:,i_channel)));

    % compute 'rfft_psd' on each frequency timeseries
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
modulation_spectrogram_data.win_size_samples = win_size;
modulation_spectrogram_data.win_shift_samples = win_shift;
modulation_spectrogram_data.n_fft_x = n_fft_x;
modulation_spectrogram_data.n_fft_y = n_fft_y;
modulation_spectrogram_data.win_funct_y = win_funct_y;
modulation_spectrogram_data.win_funct_x = win_funct_x;
modulation_spectrogram_data.n_windows = n_windows;
modulation_spectrogram_data.n_samples = spectrogram_data.n_samples;
modulation_spectrogram_data.spectrogram_data = spectrogram_data;
modulation_spectrogram_data.channel_names = channel_names;
end
