function [wcoef,wfam] = cmorlet_wavelet(x, fs, freqs, n, normalization)
%[wcoef,wfam] = cmorlet_wavelet(x, fs, freqs, n, normalization)
%CMORLET_WAVELET Wavelet transform of ONE signal
%   Continuous wavelet tranform using the Complex Morlet Wavelet
%
%   WCOEF has the same class as X
% 
%   WCOEF = CMORLET_WAVELET(SIGNAL, FS, FREQS, N, NORMALIZATION) 
%   obtains the Wavelet transform using the Complex Morlet wavelet 
%   as mother wavelet using the definition 
%   cmor = A * exp(1i*2*pi*f*time)* exp(-time.^2./(2*s^2))
%   with s = n / (2*pi*f) 
%   
%   if NORMALIZATION  is absent or equal to TRUE:
%   A = 1 / ( sqrt(szx)*((s(iwav)^2)*pi)^(1/4) );
% 
%   X     = real-valued signal or set of signals [n_samples, n_channels]
%   FS    = Sampling frequency
%   FREQS  = Array with the frequencies where the wavelet convolution is computed
%   n     = Relation between the number of cicles allowed by the Gaussian
%           defaul n = 6 
%   NORMALIZATION = Defines the value of A
%                   default normalization = true
%   WCOEF =  complex coefficients of the wavelet transform for frequencies in FREQS 
%            WCOEF = [n_samples, numel(FREQS),n_channels]
%   WFAM  =  wavelet family (set of wavlets) used. Size [n_samples, numel(FREQS)]
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
%   Raymundo Cassani.
%   June 2014

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

% Number of samples
n_samples = size(x , 1);
% Number of channels
n_channels = size(x , 2);
% Number of Wavelets
n_freqs = numel(freqs);

% Number of samples for Wavetet family
% This is equal to the number of samples needed to represent 2n cycles 
% of a sine with frequency = fres(1)[Hz], sampled at fs [Hz]. 
% This is done to ensure that every wavelet in the wavalet family will be 
% close to 0 in the negative and positive edges
n_samples_wav = round( (2 * n / freqs(1)) * fs);

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
    s(iwav) = n/(2*pi*freqs(iwav));
    gaussian_win = exp(-time.^2./(2*s(iwav)^2));   
    sinwave = exp(2*pi*1i*freqs(iwav).*time); 
    if normalization
        % each wavelet has unit energy sum(abs(wavelet).^2)) = 1
        A = 1 / ( sqrt(fs)*((s(iwav)^2)*pi)^(1/4) );
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
