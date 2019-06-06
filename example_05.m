%% Example 04
% This example shows the use of the transfomrs and their inverses 
%
% rfft()                            Fourier transform for real-valued signals 
% irfft()                           Inverse Fourier transform for real-valued signals 
%
% rfft_psd()                        Computes PSD data from x(f)
% irfft_psd()                       Recovers x(t) from its PSD data
%
% strfft_spectrogram()              Compute Spectrogram data using STFFT
% istrfft_modulation_spectrogram()  Recovers x(t) from its STFFT Spectrogram data 
%
% wavelet_spectrogram()             Compute Spectrogram using CWT
% iwavelet_modulation_spectrogram() Recovers x(t) from its CWT Spectrogram data
%

%% signal 
fs = 2000;
time_v = (0:(10*fs)-1)'/fs;
f1 = sin(2 * pi * 100 * time_v);
f2 = sin(2 * pi * 325 * time_v);
m1 = sin(2 * pi * 5 * time_v) + 2;
m2 = sin(2 * pi * 3 * time_v) + 2;

xi = [f1.*m1 + f2.*m2];
n = size(xi,1);

figure()
plot_signal(xi, fs)

%% time <--> frequency 
%% FFT and IFFT of a real-valued signal
xi_rfft = rfft(xi);
xo = irfft(xi_rfft, n);

figure()
fi = subplot(2,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(2,1,2);
plot_signal(xo, fs, 'Recovered x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );

%% PSD data obtained with rFFT, and its inverse
xi_psd = rfft_psd(xi, fs);
xo = irfft_psd(xi_psd);

figure()
fi = subplot(4,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(4,1,4);
plot_signal(xo, fs, 'Recovered x(t)')
subplot(4,1,[2,3]);
plot_psd_data(xi_psd)
title('PSD of x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );


%% time <--> time-frequency 
%% STFT Spectrogram
xi_strfft = strfft_spectrogram(xi, fs, round(fs * 0.1), round(fs * 0.05));
xo = istrfft_spectrogram(xi_strfft);

figure()
fi = subplot(4,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(4,1,4);
plot_signal(xo, fs, 'Recovered x(t)')
subplot(4,1,[2,3]);
plot_spectrogram_data(xi_strfft)
title('STFT Spectrogram of x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );

%% CWT Spectrogram
xi_cwt = wavelet_spectrogram(xi, fs);
xo = iwavelet_spectrogram(xi_cwt);

figure()
fi = subplot(4,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(4,1,4);
plot_signal(xo, fs, 'Recovered x(t)')
subplot(4,1,[2,3]);
plot_spectrogram_data(xi_cwt)
title('STFT Spectrogram of x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );


%% time <--> frequency-modulation-frequency 
%% STFT Modulation Spectrogram
xi_mod_strfft = strfft_modulation_spectrogram(xi, fs, round(fs * 0.1), round(fs * 0.05));
xo = istrfft_modulation_spectrogram(xi_mod_strfft);

figure()
fi = subplot(4,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(4,1,4);
plot_signal(xo, fs, 'Recovered x(t)')
subplot(4,1,[2,3]);
plot_modulation_spectrogram_data(xi_mod_strfft, [], [0, 1000], [0, 10])
title('STFT Modulation Spectrogram of x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );

%% cwt Modulation Spectrogram
xi_mod_cwt = wavelet_modulation_spectrogram(xi, fs);
xo = iwavelet_modulation_spectrogram(xi_mod_cwt);

figure()
fi = subplot(4,1,1);
plot_signal(xi, fs, 'Original x(t)')
fo = subplot(4,1,4);
plot_signal(xo, fs, 'Recovered x(t)')
subplot(4,1,[2,3]);
plot_modulation_spectrogram_data(xi_mod_cwt, [], [0, 1000], [0, 10])
title('STFT Modulation Spectrogram of x(t)')
linkaxes([fi, fo], 'xy')
r = corrcoef(xi, xo);
fprintf('Correlation: %0.3f \r\n', r(2) );
