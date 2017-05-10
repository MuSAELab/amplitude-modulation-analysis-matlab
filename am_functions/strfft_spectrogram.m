function spectrogram_struct = strfft_spectrogram(x, fs, win_size, win_shift, n_fft, win_funct, channel_names)
%spectrogram_struct = STRFFT_SPECTROGRAM(x, fs, win_size, win_shift, n_fft, win_funct, channel_names)
% This function computes the Short Time real FFT Spectrogram for one or a set of REAL signals 'x'
% rFFT_spectrogram and pwr_spectrogram have the same numeric class as 'x' (single, double, etc.)
%
% INPUTS:
%  x             Real-valued signal or set of signals [n_samples, n_channels]
%  fs            Sampling frequency (Hz)
%  win_size      Size of the sliding window for STFFF (samples)
%  win_shift     Shift between consecutive windows (samples)
% Optional:
%  n_fft         Number of elements to perform STFFT (default: win_size)
%  win_funct     Window to apply to the data in 'x' (default: 'blackmanharris' )
% channel_names  Cell array with names of channels (default: {'CH1', 'CH2', ... , 'CHn_channels'})
%
% OUTPUTS
%  spectrogram_struct. Output structure
%   rFFT_spectrogram  rFFT values for each window (u)
%   pwr_spectrogram   PSD values for each window (u^2 / Hz)
%   fs                Sampling frequency (Hz);
%   freq_axis         Frequency axis for rFFT_ and pwr_spectrogram (Hz)
%   freq_delta        Frequency axis step (Hz)
%   time_axis         Time axis for rFFT_ and pwr_spectrogram (s)
%   time_delta        Time axis step (s)
%   win_size          Size of the slide window (samples)
%   win_shift         Shift between consecutive windows (samples)
%   win_funct         Window function applied to each ST window the data
%   n_fft             Number of elements utilized to perform STFFT
%   n_samples_win     Number of samples in each ST window
%   n_windows         Number of ST windows
%   n_samples         Number of samples of the signal or signals 'x'
%   channel_names     cell array with names of channels
%
% Raymundo Cassani
% Nov 2016

% get class of x
x_class = class(x);

% validate 'window_funct' argument
if ~exist('win_funct','var') || isempty(win_funct)
    win_funct = 'blackmanharris';
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
        channel_names{ic} = sprintf('CH-%02d',ic);
    end
end

% compute PSD per window
for i_window = 1 : n_windows
    % ith epoch of the signal or signals
    x_epoch = squeeze(x_epoched(:, :, i_window));
    psd_struct = rfft_psd(x_epoch, fs, n_fft, win_funct, channel_names);

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

% output 'spectrogram_struct' structure
spectrogram_struct.rFFT_spectrogram = rFFT_spectrogram;
spectrogram_struct.pwr_spectrogram = pwr_spectrogram;
spectrogram_struct.fs = fs;
spectrogram_struct.freq_axis = f_ax;
spectrogram_struct.freq_delta = f_delta;
spectrogram_struct.time_axis = t_ax;
spectrogram_struct.time_delta = t_delta;
spectrogram_struct.win_size  = win_size;
spectrogram_struct.win_shift = win_shift;
spectrogram_struct.n_fft = n_fft;
spectrogram_struct.n_samples_win = n_samples_win;
spectrogram_struct.win_funct = win_funct;
spectrogram_struct.n_windows = n_windows;
spectrogram_struct.n_samples = n_samples;
spectrogram_struct.channel_names = channel_names;

end
