%% Example 03
% Script to show the computation of several Modulation Spectrograms from a signal
%
% The Modulation Spectrogram are computed with the Wavelet transform
% 
% Method A.
%   The signal, is segmented an a Modulation Spectrogram computed per Segment
%
% Method B.
%   The Spectrogram of the Full signal is obtained, after this Spectrogram
%   is segmented, and from each Segment the Modulation Spectrogram is
%   derived
%
%
% Moreover, this script compares diverse ways to compute the power of a
% signal in the Time, Time-Frequency and Frequency-Frequency domains
%
close all
clear
clc;

%% ECG data (1 channel) 
load('.\data\ecg_data.mat');

%% Segment parameters
segment_length  = 5; % seconds
segment_overlap = 0; % seconds

%% Power in Time Domain
T = 1/fs;
n = numel(x);

% Energy of the signal
energy_x = T .* sum(x.^2);
duration = T .* n;

% Power of the signal
power_x = energy_x ./ duration;

%% A. Epoching > Spectrogram > Modulation
tic

% Epochiong data 
x_segmented = epoching(x, segment_length * fs );
n_segments  = size(x_segmented, 3);

% For each Segment, compute its Spectrogram and Modulation Spectrogram 
for i_segment = 1 : n_segments
    x_tmp_wavelet = squeeze(x_segmented(:,:, i_segment));
    
    % Wavelet-based Spectrogram and Modulation Spectrogram 
    wavelet_spect_struct_a(i_segment) = wavelet_spectrogram(x_tmp_wavelet, fs);
    wavelet_modul_struct_a(i_segment) = wavelet_modspectrogram(x_tmp_wavelet, fs);   
end
toc

%% B. Spectrogram > Epoching > Modulation
tic
% Spectrogram of the Full Signal with STFFT and Wavelets
wavelet_spect_struct = wavelet_spectrogram(x ,fs);

% Epoching the Spectrogram
wavelet_spect_segmented = epoching(wavelet_spect_struct.pwr_spectrogram, segment_length * fs);
n_segments = size(wavelet_spect_segmented, 3);

% The Spectograms are scaled to represent the power of the full signal
wavelet_spect_segmented  = n_segments * wavelet_spect_segmented;

% Initialize structures, for sake of plotting
wavelet_spect_struct_b  = wavelet_spect_struct_a;
wavelet_modul_struct_b = wavelet_modul_struct_a;

% From each Segment of the Spectrogram, compute the Modulation Spectrogram
for i_segment = 1 : n_segments
    wavelet_spect_struct_b(i_segment).pwr_spectrogram = wavelet_spect_segmented(:,:,i_segment);

    % Square Root is obtained to work with the Instantaneous Amplitude
    x_tmp_wavelet = sqrt(squeeze(wavelet_spect_segmented(:,:, i_segment )));
    
    % PSD of the Spectrogram Segment
    mod_psd_wavelet = rfft_psd(x_tmp_wavelet, fs);
    
    % Place results in corresponding structure
    wavelet_modul_struct_b(i_segment).pwr_modspec = mod_psd_wavelet.PSD' / mod_psd_wavelet.freq_delta; 
end 

toc
%% Comparison
% One segment is randomly chosen 
random_segment = randi(n_segments);

for i_segment = 1 : n_segments
    if i_segment == random_segment
        figure()
        subplot(1,2,1)
        plot_spectrogram_struct(wavelet_spect_struct_a(i_segment), 1);
        subplot(1,2,2)
        plot_spectrogram_struct(wavelet_spect_struct_b(i_segment), 1);

        figure()
        subplot(1,2,1)
        plot_modspectrogram_struct(wavelet_modul_struct_a(i_segment), 1);
        subplot(1,2,2)
        plot_modspectrogram_struct(wavelet_modul_struct_b(i_segment), 1); 
    end

    pwr_spect_wavelet_a(i_segment) = sum(sum(wavelet_spect_struct_a(i_segment).pwr_spectrogram) * wavelet_spect_struct_a(1).freq_delta * wavelet_spect_struct_a(1).time_delta) ;
    pwr_spect_wavelet_b(i_segment) = sum(sum(wavelet_spect_struct_b(i_segment).pwr_spectrogram) * wavelet_spect_struct_b(1).freq_delta * wavelet_spect_struct_b(1).time_delta) ;

    pwr_modul_wavelet_a(i_segment) = sum(sum(wavelet_modul_struct_a(i_segment).pwr_modspec) * wavelet_modul_struct_a(1).freq_delta * wavelet_modul_struct_a(1).modfreq_delta);
    pwr_modul_wavelet_b(i_segment) = sum(sum(wavelet_modul_struct_b(i_segment).pwr_modspec) * wavelet_modul_struct_b(1).freq_delta * wavelet_modul_struct_b(1).modfreq_delta);
end

%% Power comparison Spectrogram and Modulation Spectrogram
figure()
title('Total Power per Epoch, based on Spectrogram')
plot([pwr_spect_wavelet_a', pwr_spect_wavelet_b'])
legend('Wavelet Spectrogram A', 'Wavelet Spectrogram B')
mean([pwr_spect_wavelet_a', pwr_spect_wavelet_b'])

figure()
title('Total Power per Epoch, based on Modulation Spectrogram')
plot([pwr_modul_wavelet_a', pwr_modul_wavelet_b'])
legend('Wavelet Modulation Spectrogram A', 'Wavelet Modulation Spectrogram B')
mean([pwr_modul_wavelet_a', pwr_modul_wavelet_b'])
