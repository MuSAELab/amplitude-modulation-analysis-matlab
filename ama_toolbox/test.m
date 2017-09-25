clc;
clear all
close all
fs = 256;
t_5s = 0 : 1/fs: 20- (1/fs);
x = [sin(2*pi*8*t_5s'),sin(2*pi*25*t_5s') ];

tic
w = wavelet_modulation_spectrogram(x, fs);
toc
tic
f = strfft_modulation_spectrogram(x, fs, 1*fs, 0.5*fs);
toc

plot_modulation_spectrogram_data(w)
plot_spectrogram_data(w.spectrogram_data)

plot_modulation_spectrogram_data(f)
plot_spectrogram_data(f.spectrogram_data)


%%
%
load('data_N030.mat');
x = eeg(1:fs*1000, 1:5);
tic
w = wavelet_modulation_spectrogram(x, fs);
toc
tic
f = strfft_modulation_spectrogram(x, fs, 1*fs, 0.5*fs);
toc


% %% 
% clc;
% clear
% n = logspace(1,5,20);
% 
% for i_n = 1 : numel(n)
%     a = [1:n(i_n)]';
%     tic
%     c = conv(a, a);
%     tm(i_n) = toc;
%     tic
%     c = conv_fft(a, a);
%     tr(i_n) = toc;
% end
% 
% loglog(n, [tm', tr'])
% legend('MATLAB', 'Ray')
% 
% %%
% %% 
% clc;
% clear
% n = logspace(1,6,20);
% for i_n = 1 : numel(n)
%     a = [1:n(i_n)];
%     tic    
%     for ix = 1 : 2
%         tc = conv(a, a);
%         if ix == 1
%         c = zeros(numel(tc),2);
%         end
%         c(:, i_n) = tc;
%     end
%     tm(i_n) = toc;
%     tic
%     c = conv_fft(a', [a', a']);
%     tr(i_n) = toc;
%     
%     tic
%     c = conv_fft2(a', [a', a']);
%     tr2(i_n) = toc;
%     
%     
%     
% end
% 
% semilogx(n, [tm', tr', tr2'])
% legend('MATLAB', 'Ray', 'Ray2')


