function spectrogram_data = strfft_spectrogram(x, fs, win_size, win_shift, n_fft, win_function, channel_names)
%spectrogram_data = STRFFT_SPECTROGRAM(x, fs, win_size, win_shift, n_fft, win_function, channel_names)
%     Compute the Short Time real FFT Spectrogram for one or a set of REAL signals 'x'.
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
%     n_fft : Number of samples to compute the FFT
%             (Default = n_samples in array x)   
%     win_function : Window function applied to the signal 
%         (Default 'Hamming')
%     channel_names : Names of the signals
%         (Default Signal-XX with XX 1, 2, ... n_channels) 
% 
%     Returns
%     -------
%     spectrogram_data : Structure with Spectrogram data, with the elements:
%        rFFT_spectrogram
%            rFFT values for each window (u), scaled by the Window RMS       
%        power_spectrogram :
%            PSD values for each window (u^2 / Hz) 
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
%        win_size_samples :
%            Size of the sliding window for STFFF (samples)
%        win_shift_samples :
%            Shift between consecutive windows (samples)   
%        n_fft :
%            Number of elements utilized to perform FFT    
%        win_function :
%            Window applied to the data in 'x'           
%        n_windows :
%            Number of ST windows
%        n_samples :
%            Number of samples of the signal or signals 'x'
%        channel_names 
%            Names of channels
           
% get class of x
x_class = class(x);

% validate 'window_funct' argument
if ~exist('win_function','var') || isempty(win_function)
    win_function = 'hamming';
end

% validate 'n_fft' argument
if ~exist('n_fft','var') || isempty(n_fft)
    n_fft = win_size;
end

% validate 'channel_names' argument
if ~exist('channel_names','var') || isempty(channel_names)
    channel_names = {};
end

% round win_size and win_shift
win_size = round(win_size);
win_shift = round(win_shift);

% time axis step for Spectrogram
t_delta = win_shift / fs;

% number of samples for original signals
n_samples = size(x , 1);

% create time vector 'time_vct' for signal 'x'
time_vct = (0 : n_samples - 1) / fs;

% epoch signal or signals 'x'
[x_epoched, ~, ix] = epoching(x, win_size, win_size - win_shift);

% time axis for Spectrogram
t_ax = time_vct(ix);

% spectrogram parameters
n_windows  = size(x_epoched, 3);
n_channels = size(x_epoched, 2);
n_samples_win  = size(x_epoched, 1);

% generate default channel names, if needed
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('Signal-%02d',ic);
    end
end

% compute PSD per window
for i_window = 1 : n_windows
    % ith epoch of the signal or signals
    x_epoch = squeeze(x_epoched(:, :, i_window));
    psd_struct = rfft_psd(x_epoch, fs, n_fft, win_function, channel_names);

    % initialize arrays for spectrogram data
    if i_window == 1
        % frequency Axis for spectrogram
        f_ax = psd_struct.freq_axis;
        % delta Frequency
        f_delta = psd_struct.freq_delta;
        % initialize 'rFFT_spectrogram' and 'pwr_spectrogram'
        rFFT_spectrogram = zeros(n_windows, numel(f_ax), n_channels, x_class);
        pwr_spectrogram  = zeros(n_windows, numel(f_ax), n_channels, x_class);
    end
    % rFFT data
    rFFT_spectrogram(i_window, :, :) = psd_struct.rFFT;
    % power data
    pwr_spectrogram(i_window, :, :) = psd_struct.PSD;
end

% scale 'pwr_spectrogram' by number of windows and time delta
pwr_spectrogram = pwr_spectrogram / (n_windows * t_delta);

% output 'spectrogram_data' structure
spectrogram_data.rFFT_spectrogram = rFFT_spectrogram;
spectrogram_data.power_spectrogram = pwr_spectrogram;
spectrogram_data.fs = fs;
spectrogram_data.freq_axis = f_ax;
spectrogram_data.freq_delta = f_delta;
spectrogram_data.time_axis = t_ax;
spectrogram_data.time_delta = t_delta;
spectrogram_data.win_size_samples  = win_size;
spectrogram_data.win_shift_samples = win_shift;
spectrogram_data.n_fft = n_fft;
spectrogram_data.win_function = win_function;
spectrogram_data.n_windows = n_windows;
spectrogram_data.n_samples = n_samples;
spectrogram_data.channel_names = channel_names;

end
