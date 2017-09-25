function spectrogram_data = wavelet_spectrogram(x, fs, n_cycles, freq_vct, channel_names)
%spectrogram_struct = WAVELET_SPECTROGRAM(x, fs, n_cycles, freq_vct, channel_names)
%     Compute the Spectrogram using the Complex Morlet wavelet for one or a set of REAL signals 'x'. 
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
%     channel_names : Names of the signals
%         (Default Signal-XX with XX 1, 2, ... n_channels) 
% 
%     Returns
%     -------
%     spectrogram_data : Structure with Spectrogram data, with the elements:
%        wavelet_coefficients
%            Coefficients of the Wavelet transformation (u)       
%        power_spectrogram :
%            Power spectrogram (u^2 / Hz) 
%        fs : 
%            Sampling frequency (Hz)
%        freq_axis :
%            Frequency axis for rFFT and PSD (Hz)
%        freq_delta :
%            Frequency axis step (Hz)
%        time_axis :
%            Time axis for rFFT_spectrogram and power_spectrogram (s)       
%        time_delta :
%            Time axis step (s)
%        n_cycles : 
%            Number of cicles used inside the Gaussian curve 
%        wavelet_kernels :
%            Wavelet kernels used to obtain the wavelet coefficients
%        n_samples :
%            Number of samples of the signal or signals 'x'
%        channel_names :
%            Names of channels

% validate 'freq_vct' argument
if ~exist('freq_vct','var') || isempty(freq_vct)
    freq_vct = 1 : floor(fs / 2);
end

% validate 'n_cycles' argument
if ~exist('n_cycles','var') || isempty(n_cycles)
    n_cycles = 6;
end

% validate 'channel_names' argument
if ~exist('channel_names','var') || isempty(channel_names)
    channel_names = {};
end

% time delta
t_delta = 1 / fs;
% frequency delta
f_delta = freq_vct(2) - freq_vct(1);

% create time vector 'time_vct' for signal 'x'
time_vct = (0 : size(x,1) - 1) / fs;

% time axis for spectrogram
t_ax = time_vct;
% frequency axis for spectrogram
f_ax = freq_vct;

% number of channels and number of samples
n_channels = size(x, 2);
n_samples  = size(x, 1);

% generate default channel names
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('CH-%02d',ic);
    end
end

% wavelet transform
[wcoef,wfam] = cmorlet_wavelet(x, fs, freq_vct, n_cycles, true);
% power from Wavelet coefficients
pwr_spectrogram = abs(wcoef).^2;
pwr_spectrogram = pwr_spectrogram * 2 / (fs * n_samples);

% output 'spectrogram_data' structure
spectrogram_data.wavelet_coefficients = wcoef;
spectrogram_data.power_spectrogram = pwr_spectrogram;
spectrogram_data.fs = fs;
spectrogram_data.freq_axis = f_ax;
spectrogram_data.freq_delta = f_delta;
spectrogram_data.time_axis = t_ax;
spectrogram_data.time_delta = t_delta;
spectrogram_data.n_cycles = n_cycles;
spectrogram_data.wavelet_kernels = wfam;
spectrogram_data.n_samples = n_samples;
spectrogram_data.channel_names = channel_names;

end
