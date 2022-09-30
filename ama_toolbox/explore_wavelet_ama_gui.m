function explore_wavelet_ama_gui(X, fs, Names, c_map)
% Analysis of a Signal in Frequency-Frequency Domain
% Time -> Time-Frequency transformation performed with Wavelet Transform (Complex Morlet)
%
% INPUTS:
%  X     Real-valued column-vector signal or set of signals [n_samples, n_channels]
%  fs    Sampling frequency (Hz)
% Optional:
%  Names (Optional) Name of the signal(s), String or Cell Array (Strings)
%  c_map (Optional) Colormap, Default 'viridis'
%

clc;

% Global variable with the value of the Key Pressed
global key_pressed
key_pressed = '';

%% verify input arguments
if ~exist('X','var')
    error('Variable X is requiered');
end

if ~exist('fs','var')
    error('Variable fs is requiered');
end

% number of channels
n_channels = size(X, 2);

% Validate 'c_map' argumet
if ~exist('c_map','var') || isempty(c_map)
    c_map = 'viridis';
end

% verify channel names
if exist('Names','var')
    if ischar(Names) && n_channels == 1
        Names = {Names};
    end
end
if ~exist('Names','var') || numel(Names) ~= n_channels
    Names = cell(1, n_channels);
    for ix = 1 : n_channels
        Names{ix} = sprintf('Signal-%02d', ix);
    end
end

% check if 'am_functions' folder exist, if does, add it to the MATLAB path
if exist('./am_functions', 'dir')
    addpath('./am_functions');
end

% verify dependencies
if ~exist('wavelet_modulation_spectrogram.m','file')
    error('Dependencies for explore_wavelet_am_gui were not satisfied');
end

%% Amplitude Modulation Analysis
% default Modulation Analysis parameters
n_cycles     =  6;    % number of cycles (for Complex Morlet)
seg_size_sec =  8;    % segment of signal to compute the Modulation Spectrogram (seconds)
seg_shft_sec =  8;    % shift between consecutive segments (seconds)
freq_range   = [];    % limits [min, max] for the conventional frequency axis (Hz)
mfreq_range  = [];    % limits [min, max] for the modulation frequency axis (Hz)
freq_color   = [];    % limits [min, max] for the power in Spectrogram (dB)
mfreq_color  = [];    % limits [min, max] for the power in Modulation Spectrogram (dB)

% initial channel and segment
ix_channel = 1;
ix_segment = 1;

% other variables
n_segments = [];
h_ts = [];
h_tf = [];
h_area1 = [];
h_area2 = [];
x_segments = [];

%% Live GUI
h_fig = figure('Name', 'Explore Wavelet Amplitude Modulation', 'KeyPressFcn', {@key_pressed_fcn});
first_run();

while true
   pause(0.01);
   % handles user request in GUI
    switch key_pressed
        case 'escape' %ESC: Exit
            break
        case {'leftarrow', 'left'} %Left arrow: Previous Segment
            ix_segment = ix_segment - 1;
            key_pressed = '';
            update_plots();
        case {'rightarrow', 'right'} %Right arrow: Next Segment
            ix_segment = ix_segment + 1;
            key_pressed = '';
            update_plots();
        case {'uparrow', 'up'} %Up arrow: Previous
            ix_channel = ix_channel - 1;
            key_pressed = '';
            first_run();
        case {'downarrow', 'down'} %Down arrow: Next Channel
            ix_channel = ix_channel + 1;
            key_pressed = '';
            first_run();
        case {'a', 'A'}    %A: Back 5 Segments
            ix_segment = ix_segment - 5;
            key_pressed = '';
            update_plots();
        case {'d', 'D'}    %D: Advance 5 Segment
            ix_segment = ix_segment + 5;
            key_pressed = '';
            update_plots();
        case {'w', 'W'}    %W  Previous 5 channels
            ix_channel = ix_channel - 5;
            key_pressed = '';
            first_run();
        case {'s', 'S'}    %S  Next 5 Channels
            ix_channel = ix_channel + 5;
            key_pressed = '';
            first_run();
        case {'u', 'U'}
            key_pressed = '';
            update_parameters();
    end
end

    % executed at first time running or when an update in parameters is performed
    function first_run()
        ix_segment = 1;
        clc;
        disp('computing full-signal spectrogram...');

        % constrain ix_channel to [1 : n_channels]
        ix_channel = max([1, ix_channel]);
        ix_channel = min([n_channels, ix_channel]);

        % wavelet modulation spectrogram parameters in samples
        seg_size_smp = round(seg_size_sec * fs);  % (samples)
        seg_shft_smp = round(seg_shft_sec * fs);  % (samples)

        % signal for analysis
        x_probe = X(:, ix_channel);

        % segment of signal under analysis
        x_segments = epoching(x_probe, seg_size_smp, seg_size_smp - seg_shft_smp);
        n_segments = size(x_segments,3);

        % plot complete time series
        h_ts = subplot(4,2,[1,2]);
        plot_signal(x_probe(:), fs, Names{ix_channel});
        colorbar();
        time_lim = get(gca, 'XLim' );
        h_area1 = [];

        % compute and plot complete spectrogram
        x_spectrogram = wavelet_spectrogram(x_probe, fs, n_cycles, [], Names(ix_channel));
        h_tf = subplot(4,2,[3,4]);
        plot_spectrogram_data(x_spectrogram, [], [], freq_range, freq_color, c_map)
        set(gca, 'XLim', time_lim);
        h_area2 = [];

        % update plots
        update_plots();
    end

    % update plots when a new segment is selected
    function update_plots()
        % constrain ix_segment to [1 : n_segments]
        ix_segment = max([1, ix_segment]);
        ix_segment = min([n_segments, ix_segment]);

        % delete old areas (if existent)
        if ~isempty(h_area1)
            delete(h_area1);
            delete(h_area2);
        end

        % select segment
        x = x_segments(:, :, ix_segment);

        % compute and plot Modulation Spectrogram
        clc;
        disp('computing modulation spectrogram...');
        x_wavelet_modspec = wavelet_modulation_spectrogram(x, fs, n_cycles, [], 2, [], Names(ix_channel));
        subplot(4,2,[6,8])
        plot_modulation_spectrogram_data(x_wavelet_modspec, [], freq_range, mfreq_range, mfreq_color,c_map)
        % Uncomment for log axes
        %set(gca,'XScale','log');
        %set(gca,'YScale','log');

        % plot time series for segment
        subplot(4,2,5)
        plot_signal(x, fs, Names{ix_channel});
        colorbar();
        time_lim = get(gca, 'XLim' );

        % plot spectrogram for segment
        subplot(4,2,7)
        plot_spectrogram_data(x_wavelet_modspec.spectrogram_data, [], [], freq_range, freq_color,c_map)
        set(gca, 'XLim', time_lim);

        % highlight area under analysis in time series
        seg_ini_sec = (ix_segment - 1) * seg_shft_sec;
        seg_end_sec = seg_ini_sec + seg_size_sec;
        subplot(h_ts)
        h_area1 = varea([seg_ini_sec, seg_end_sec ],{'g'});
        % highlight area under analysis in Spectrogram
        subplot(h_tf)
        h_area2 = varea([seg_ini_sec, seg_end_sec ],{'g'}, 0.8);
        disp('done!');

        % display information about analysis
        clc;
        fprintf('signal name            : %s\n', Names{ix_channel} );
        fprintf('segment size  (seconds): %0.3f\n', seg_size_sec);
        fprintf('segment shift (seconds): %0.3f\n', seg_shft_sec);
        fprintf('segment position  (sec): %0.3f\n', seg_ini_sec)
        fprintf('n cycles Complex Morlet: %d\n', n_cycles);

        drawnow();
    end

    % update parameters for AM analysis
    function update_parameters()
        % GUI to get new parameters
        prompt = {'Segment         (seconds): ', ...
                  'Segment shift   (seconds): ', ...
                  'N Cycles                 : ', ...
                  'Freq Conv. min,Max  (Hz) : ', ...
                  'Spectr Pwr min,Max  (dB) : ', ...
                  'Freq Mod.  min,Max  (Hz) : ', ...
                  'ModSpec Pwr min,Max (dB) : ', ...
                  };

        win_name = 'Modulation Analysis Parameters'; numlines = 1;

        defaultanswer = {num2str(seg_size_sec), num2str(seg_shft_sec), ...
                         num2str(n_cycles), ...
                         num2str(freq_range), num2str(freq_color), ...
                         num2str(mfreq_range), num2str(mfreq_color),...
                         };

        options.Resize = 'off';

        answer = inputdlg(prompt,win_name,numlines,defaultanswer,options);

        % if the user clicks the 'Cancel' button, 'answer' is empty
        if ~isempty(answer)
            seg_size_sec = str2double(answer{1});  % (seconds)
            seg_shft_sec = str2double(answer{2});  % (seconds)
            n_cycles     = str2double(answer{3});  % (cycles)
            freq_range   = validate_range(answer{4});  % (Hz)
            freq_color   = validate_range(answer{5});  % (dB)
            mfreq_range  = validate_range(answer{6});  % (Hz)
            mfreq_color  = validate_range(answer{7});  % (dB)

            % call first_run()
            first_run();
        end
    end

    function out = validate_range(inStr)
        out = [];
        strParts = strsplit(inStr);
        if length(strParts) ~= 2
            return
        end
        numParts = str2double(strParts);
        if any(isnan(numParts))
            return
        end
        out = numParts;
    end

% if the user press ESC, close the GUI and clean console
close(h_fig);
clc;

end
