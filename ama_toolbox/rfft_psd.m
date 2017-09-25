function psd_data = rfft_psd(x, fs, n_fft, win_funct, channel_names)
%psd_struct = RFFT_PSD(x, fs, n_fft, win_funct, channel_names)
%
%     Compute the PSD for one or a set of REAL signals.
%         
%     Parameters
%     ----------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_channels)
%     fs : Sampling frequency 
%         in Hz
%     n_fft : Number of samples to compute the FFT
%             (Default = n_samples in array x)   
%     win_funct : Window function applied to the signal 
%         (Default 'Hamming')
%     channel_names : Names of the signals
%         (Default Signal-XX with XX 1, 2, ... n_channels) 
% 
%     Returns
%     -------
%     psd_data : Structure with PSD data, with the elements:
%        rFFT
%            First half of the FFT(x) (u), scaled by the Window RMS       
%        PSD
%            Power Spectrum Density (u^2 / Hz) 
%        fs
%            Sampling frequency (Hz)
%        freq_axis
%            Frequency axis for rFFT and PSD (Hz)
%        freq_delta
%            Frequency axis step (Hz)
%        n_samples
%            Number of samples of the signal or signals 'x'
%        n_fft
%            Number of elements utilized to perform FFT
%        win_funct
%            Window applied to the data in 'x'
%        channel_names 
%            Names of channels


% validate 'n_fft' argument
if ~exist('n_fft','var') || isempty(n_fft)
    n_fft = size(x,1);
end

% validate 'win_name' argument
if ~exist('win_funct','var') || isempty(win_funct)
    win_funct = 'hamming';
end

% validate 'channel_names' argument
if ~exist('channel_names','var') || isempty(channel_names)
    channel_names = {};
end

% number of channels and number of samples
[n_samples, n_channels] = size(x);

% generate default channel names, if needed
if isempty(channel_names)
    for ic = 1 : n_channels
        channel_names{ic} = sprintf('Signal-%02d',ic);
    end
end

% windowing data
win = window(win_funct, n_samples);
win_rms = sqrt(sum(win.^2) / n_samples);
win_mat = repmat(win, 1, n_channels);
x = x.*win_mat;

% real FFT with zero padding if n_fft ~= n_samples
X = rfft(x, n_fft);
% spectrum scaled by window RMS
X = X / win_rms;
% power spectrum
X_pwr = X .* conj(X);
X_pwr = X_pwr * (1/n_fft.^2);

% adjust for even and odd number of elements
if mod(n_fft,2)
    % odd case
    n_freqs = (n_fft + 1) / 2;
    % double all frequency components except DC component
    X_pwr(2:end,:) = X_pwr(2:end,:) * 2;
else
    % even case
    n_freqs = (n_fft / 2) + 1;
    % double all frequency components except DC and fs/2 components
    X_pwr(2:end-1,:) = X_pwr(2:end-1,:) * 2;
end

% frequency axis step
f_delta = (fs / n_fft);

% scale PSD with the frequency step
psd = X_pwr ./ f_delta;

% frequency axis for spectrum
f_ax = (0: n_freqs-1) * f_delta;

% power of the signal can be computed as:
% power_x = sum(psd) * f_delta == sum(x.^2) / n_samples
% which agrees with Parseval's Theorem

% output 'psd_struct' structure
psd_data.rFFT = X;
psd_data.PSD = psd;
psd_data.fs = fs;
psd_data.freq_axis = f_ax;
psd_data.freq_delta = f_delta;
psd_data.n_samples = n_samples;
psd_data.n_fft = n_fft;
psd_data.win_funct = win_funct;
psd_data.channel_names = channel_names;
