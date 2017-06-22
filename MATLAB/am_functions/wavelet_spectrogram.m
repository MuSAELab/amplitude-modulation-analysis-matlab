function spectrogram_struct = wavelet_spectrogram(x, fs, n_cycles, freq_vct, channel_names)
%spectrogram_struct = WAVELET_SPECTROGRAM(x, fs, n_cycles, freq_vct, channel_names)
% This function computes the Spectrogram using the Complex Morlet wavelet
% for one or a set of REAL signals 'x'
% rFFT_spectrogram and pwr_spectrogram have the same numeric class as 'x'
%
% INPUTS:
%  x             Real-valued signal or set of signals [n_samples, n_channels]
%  fs            Sampling frequency (Hz)
% Optional:
%  freq_vct      Indicates the frequencies for the Wavelet transform
%                 Default = 1 : floor(fs/2);
%  n_cycles      Approx. number of cycles in the Gaussian kernel (default = 6)
% channel_names  Cell array with names of channels (default: {'CH1', 'CH2', ... , 'CHn_channels'})
%
% OUTPUTS:
%  spectrogram_struct. Output structure
%   wavelet_coef      Coefficients of the Wavelet transformation (u)
%   pwr_spectrogram   Scaled power
%   fs                Sampling frequency (Hz);
%   freq_axis         Frequency axis for rFFT_ and pwr_spectrogram (Hz)
%   freq_delta        Frequency axis step (Hz)
%   time_axis         Time axis for rFFT_ and pwr_spectrogram (s)
%   time_delta        Time axis step (s)
%   n_cycles          Number of cycles in the Gaussion kernel
%   wavelet_kernels   Wavelet kernels to obtain the wavelet coefficients
%   n_samples         Number of samples of the signal or signals 'x'
%   channel_names     Cell array with names of channels
%
% Raymundo Cassani
% April 2017

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

% output 'spectrogram_struct' structure
spectrogram_struct.wavelet_coef = wcoef;
spectrogram_struct.pwr_spectrogram = pwr_spectrogram;
spectrogram_struct.fs = fs;
spectrogram_struct.freq_axis = f_ax;
spectrogram_struct.freq_delta = f_delta;
spectrogram_struct.time_axis = t_ax;
spectrogram_struct.time_delta = t_delta;
spectrogram_struct.n_cycles = n_cycles;
spectrogram_struct.wavelet_kernels = wfam;
spectrogram_struct.n_samples = n_samples;
spectrogram_struct.channel_names = channel_names;

end
