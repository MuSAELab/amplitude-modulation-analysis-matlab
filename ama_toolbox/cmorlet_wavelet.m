function [wcoef,wfam] = cmorlet_wavelet(x, fs, freq_vct, n, normalization)
%[wcoef,wfam] = cmorlet_wavelet(x, fs, freqs, n, normalization)
%CMORLET_WAVELET Perform the continuous wavelet (CWT) tranform using the complex Morlet wavelet.
%
%     Parameters
%     ----------
%     x  : 1D array with shape (n_samples) or 2D array with shape (n_samples, n_channels)
%     fs : Sampling frequency in Hz
%     freq : 1D array with frequencies to compute the CWT
%            (Default = [1 : 1 : fs/2])
%     n : Number of cicles inside the Gaussian curve (Default 6)
%     normalization: (Default True) Scale each wavelet to have energy equal to 1
% 
% 
%     Returns
%     -------
%     wcoef : Complex wavelet coefficients 
%             2D array with shape [n_samples, n_freqs] if `x` is 1D array
%             3D array with shape [n_samples, n_freqs, n_channels] if `x` is 2D array
%     
%     wfam  : 2D array with shape [n_wavelet_samples, n_freqs] where each column
%             corresponds to the a member of the wavelet family
%    
%
% Example-1:
% Test signal [x]: 5 seconds 10Hz sine followed by 8 seconds 25Hz sine
% sampled at 256Hz
% Spectrogram is computed for frequencies from 1 to 100Hz, Normalization = true
% and n equal to 6
% 
% Code:
%  fs = 256;
%  t_5s = (0:(5*fs)-1)'/fs;
%  freqs = 1:100;
%  x = [sin(8 * 2 * pi * t_5s); sin(25 * 2 * pi * t_5s)] ;
%  t = (0:numel(x)-1)/fs;
%  wcoef = cmorlet_wavelet(x,fs,freqs);
%  surf(t,freqs,20*log10(abs(wcoef')+eps),'LineStyle','none');
%  view(0,90);
%  xlabel('seconds');
%  ylabel('Hz')
%

% Get class of x
x_class = class(x);

% Validate 'n' argument, number of cicles in the Gaussian
if ~exist('n','var') || isempty(n)
    n = 6;
end

% Validate 'normalization' argument, number of cicles in the Gaussian
if ~exist('normalization','var') || isempty(normalization)
    normalization = true;
end

% validate 'freq_vct' argument
if ~exist('freq_vct','var') || isempty(freq_vct)
    freq_vct = 1 : floor(fs / 2);
end

% Number of samples
n_samples = size(x , 1);
% Number of channels
n_channels = size(x , 2);
% Number of Wavelets
n_freqs = numel(freq_vct);

% Number of samples for Wavetet family
% This is equal to the number of samples needed to represent 2n cycles 
% of a sine with frequency = fres(1)[Hz], sampled at fs [Hz]. 
% This is done to ensure that every wavelet in the wavalet family will be 
% close to 0 in the negative and positive edges
n_samples_wav = round( (2 * n / freq_vct(1)) * fs);

% Create time vector for Wavelet family
half = floor(n_samples_wav/2);
if mod(n_samples_wav,2) == 1 % odd samples
    time = (-(half):half)/fs;
else % even samples
    time = (-(half-1):half)/fs;
end

% initialize Wavelet family matrix
wfam = zeros(numel(time), n_freqs, x_class);

% for each frequency defined in FREQ, create its respective Wavelet
for iwav = 1 : n_freqs
    s(iwav) = n/(2*pi*freq_vct(iwav));
    gaussian_win = exp(-time.^2./(2*s(iwav)^2));   
    sinwave = exp(2*pi*1i*freq_vct(iwav).*time); 
    if normalization
        % each wavelet has unit energy sum(abs(wavelet).^2)) = 1
        A = 1 / ((s(iwav)^2) * pi)^(1/4);
    else 
        A = 1;
    end
  
    % Complex morlet wavelet
    wfam(:,iwav) = A * sinwave .* gaussian_win;
end

for i_channel = 1 : n_channels
    % one channel
    x_tmp = x(: , i_channel);
    % convolution between signal X and the each Wavelt in the Wavelet family
    tmp = conv_fft(x_tmp , wfam , 'same');
    if i_channel == 1
        % initialize wcoef
        wcoef = zeros(n_samples, n_freqs, n_channels, x_class);
    end
    wcoef(: , : , i_channel) = tmp;
    
end

end
