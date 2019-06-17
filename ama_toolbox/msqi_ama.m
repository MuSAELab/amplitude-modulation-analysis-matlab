function [msqi_value, hr_value, modulation_spectrogram] = msqi_ama(x, fs)
% [msqi_value, hr_value, modulation_spectrogram] = msqi_ama(x, fs)
%
% Computes the Modulation Spectrum-Based ECG Quality Index (MSQI) for one or 
% many ECG signals defined in x, sampled with a sampling frequency fs
% 
%     Parameters
%     ----------
%     x  : 1D array with shape (n_samples) or
%          2D array with shape (n_samples, n_signals)
%     fs : Sampling frequency in Hz
% 
%     Returns
%     -------
%     msqi_value : MSQI value or values 
%     hr_value   : HR values or values
%     modulation_spectrogram : Structure or structures of modulation spectrogram
%
% The MSQI is presented in:
%
% D. P. Tobon V., T. H. Falk, and M. Maier, "MS-QI:  A  Modulation
% Spectrum-Based ECG Quality Index for Telehealth Applications", IEEE
% Transactions on Biomedical Engineering, vol. 63, no. 8, pp. 1613-1622,
% Aug. 2016

    % values for the STFFT transformation
    win_size_sec  = 0.125;   % seconds
    win_over_sec  = 0.09375; % seconds
    nfft_factor_1 = 16;
    nfft_factor_2 = 4;

    win_size_smp = round(win_size_sec * fs); % samples
    win_over_smp = round(win_over_sec * fs); % samples
    win_shft_smp = win_size_smp - win_over_smp;

    % computes modulation spectrogram
    modulation_spectrogram = strfft_modulation_spectrogram(x, fs, win_size_smp, win_shft_smp, ...
                              nfft_factor_1, 'cosinewin', nfft_factor_2, 'cosinewin' );

    % find fundamental frequency (HR)
    % f = (0, 40)Hz 
    [~, ix_f_00]  = min(abs(modulation_spectrogram.freq_axis - 0));  
    [~, ix_f_40]  = min(abs(modulation_spectrogram.freq_axis - 40));

    % look for the maximum only from 0.6 to 3 Hz (36 to 180 bpm)
    valid_f_ix = and(modulation_spectrogram.freq_mod_axis > 0.6,  modulation_spectrogram.freq_mod_axis < 3);

    % number of epochs
    n_epochs = size(modulation_spectrogram.power_modulation_spectrogram, 3);

    msqi_value = zeros(n_epochs, 1);
    hr_value   = zeros(n_epochs, 1);

    for ix_epoch = 1 : n_epochs

        B = sqrt(modulation_spectrogram.power_modulation_spectrogram(:, :, ix_epoch));

        % scale to maximun of B
        B = B / max(B(:));

        % add B in the conventional frequency axis from 0 to 40 Hz
        tmp = sum(B(ix_f_00:ix_f_40, :));

        % look for the maximum only from 0.6 to 3 Hz (36 to 180 bpm)
        tmp(~valid_f_ix) = 0;
        [~, ix_fm_max] = max(tmp);
        freq_funda = modulation_spectrogram.freq_mod_axis(ix_fm_max);

        % 'energy' in lobes
        eme = 0;
        for ix_harm = 1 : 4
            [~, ix_fm] = min(abs(modulation_spectrogram.freq_mod_axis - (ix_harm * freq_funda) ));  
            ix_b = round(.3125 / modulation_spectrogram.freq_mod_delta );   %.3125Hz, half lobe used in Diana's code
            % EME
            eme = eme + sum(sum(B(ix_f_00 : ix_f_40 , ix_fm-ix_b : ix_fm+ix_b )));
        end

        % total 'energy'
        tme = sum(B(:));

        % RME
        rme = tme - eme;

        % MS-QI
        msqi_value(ix_epoch) = eme / rme;
        hr_value(ix_epoch) = freq_funda * 60;
    end

end

