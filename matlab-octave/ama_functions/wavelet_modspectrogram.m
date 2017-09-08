function modspectrogram_struct = wavelet_modspectrogram( x, fs, n_cycles, freq_vct, fft_factor_x, win_funct_x, channel_names)
%modspect_struct = WAVELET_MODSPECTROGRAM( x, fs, n_cycles, freq_vct, fft_factor_x, win_funct_x, channel_names)
% This function computes the Modulation Spectrogram (modspectrogram) for one or a set of REAL signals 'x'
% rFFT_modspec, pwd_modspec and the spectrogram_structure have the same numeric class as 'x' (single, double, etc.)
%
% INPUTS:
%  x             Real-valued column-vector signal or set of signals [n_samples, n_channels]
%  fs            Sampling frequency (Hz)
% Optional:
%  freq_vct      Indicates the frequencies for the Wavelet transform
%                 Default = 1 : floor(fs/2);
%  n_cycles      Approx. number of cycles in the Gaussian kernel (default = 6)
%  fft_factor_x  The number of elements to perform the FFT on the spectrogram, it is given as:
%                  n_fft_x  = fft_factor_x * n_windows
%                  (default, fft_factor_x = 1)
%  win_funct_x   Window to apply in the FFT of spectrogram (default: 'blackmanharris' )
%  channel_names  Cell array with names of channels (default: {'CH1', 'CH2', ... , 'CHn_channels'})
%
% OUTPUTS:
%  modspectrogram_struct. Output structure
%   rFFT_modspec     rFFT values for the Modulation Spectrogram
%   pwd_modspec      PSD values for the Modulation Spectrogram
%   fs               Sampling frequency (Hz);
%   fs_mod           Sampling frequency of modulation-frequency (Hz)
%   freq_axis        Conventional-frequency axis for rFFT_modspec and pwr_modspec (Hz)
%   freq_delta       Conventional-frequency step (Hz)
%   modfreq_axis     Modulation-frequency axis for rFFT_modspec and pwr_modspec (Hz)
%   modfreq_delta    Modulation-frequency step (Hz)
%   n_fft_x          The number of elements to perform the FFT on the spectrogram
%   win_funct_x      Window to apply in the FFT of spectrogram (default: 'blackmanharris' )
%   n_samples        Number of samples of the signal or signals 'x'
%   channel_names    Cell array with names of channels
%  spectrogram_structure  Spectrogram structure (see wavelet_spectrogram.m)
%
%
% Raymundo Cassani
% April 2017

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
    win_funct_x = 'blackmanharris';
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
[n_windows, n_freqs, n_channels] =  size(spectrogram_data.wavelet_coef);

% number of elements for FFT of the spectrogram
n_fft_x =  fft_factor_x * n_windows;

% generate default channel names
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('CH-%02d',ic);
    end
end

% modulation sampling frequency
fs_mod = fs;

% the AM analysis is made in the Amplitude derived from the Power Spectrogram
for i_channel = 1 : n_channels
  % data to generate the Modulation Spectrogram
    spectrogram_1ch = sqrt(spectrogram_data.pwr_spectrogram(:,:,i_channel)) ;
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

% output 'modspectrogram_struct' structure
modspectrogram_struct.rFFT_modspec = rFFT_modspec;
modspectrogram_struct.pwr_modspec = pwr_modspec;
modspectrogram_struct.fs = fs;
modspectrogram_struct.fs_mod = fs_mod;
modspectrogram_struct.freq_axis = spectrogram_data.freq_axis;
modspectrogram_struct.freq_delta = spectrogram_data.freq_delta;
modspectrogram_struct.modfreq_axis = fmod_ax;
modspectrogram_struct.modfreq_delta = fmod_delta;
modspectrogram_struct.n_fft_x = n_fft_x;
modspectrogram_struct.win_funct_x = win_funct_x;
modspectrogram_struct.n_samples = spectrogram_data.n_samples;
modspectrogram_struct.channel_names = channel_names;
modspectrogram_struct.spectrogram_structure = spectrogram_data;
end
