function explore_stfft_ama_gui(X, fs, Name, c_map)
% Analysis of a Signal in Frequency-Frequency Domain
% Time -> Time-Frequency transformation performed with STFFT
%
% INPUTS:
%  X     Real-valued column-vector signal or set of signals [n_samples, n_channels]
%  fs    Sampling frequency (Hz)
% Optional:
%  Name  (Optional) Name of the signal(s), String or Cell Array (Strings)
%
% Raymundo Cassani
% April 2017

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

% verify channel names
if exist('Name','var') && ischar(Name) && n_channels == 1
    Name = {Name};
end

% generate channel names if needed
if ~exist('Name','var') || isempty(Name) || numel(Name) ~= n_channels
    Name = cell(1, n_channels);
    for ix = 1 : n_channels
    Name{ix} = sprintf('Signal-%02d', ix);
    end
end

% check if 'am_functions' folder exist, if does, add it to the MATLAB path
if exist('./am_functions', 'dir')
    addpath('./am_functions');
end

% verify dependencies
if ~exist('strfft_modulation_spectrogram.m','file');
    error('Dependencies for explore_stfft_am_gui were not satisfied');
end

%% Amplitude Modulation Analysis
% default Modulation Analysis parameters
win_size_sec = 0.5;   % window length for the STFFT
win_shft_sec = 0.02;  % shift between consecutive windows (seconds)
seg_size_sec = 8;     % segment of signal to compute the Modulation Spectrogram (seconds)
seg_shft_sec = 8;     % shift between consecutive segments (seconds)
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
name = '';
win_size_smp = [];
win_shft_smp = [];

%% Live GUI
h_fig = figure('Name', 'Explore STFFT Amplitude Modulation', 'KeyPressFcn', {@key_pressed_fcn});
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
        case {'u', 'U'};
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

        % STFFT modulation spectrogram parameters in samples
        win_size_smp = round(win_size_sec * fs);  % (samples)
        win_shft_smp = round(win_shft_sec * fs);  % (samples)
        seg_size_smp = round(seg_size_sec * fs);  % (samples)
        seg_shft_smp = round(seg_shft_sec * fs);  % (samples)

        % signal for analysis
        x_probe = X(:, ix_channel);
        name = Name(ix_channel);

        % segment of signal under analysis
        x_segments = epoching(x_probe, seg_size_smp, seg_size_smp - seg_shft_smp);
        n_segments = size(x_segments,3);

        % plot complete time series
        h_ts = subplot(4,2,[1,2]);
        plot_signal(x_probe(:), fs, name);
        colorbar();
        time_lim = get(gca, 'XLim' );

        % compute and plot complete spectrogram
        x_spectrogram = strfft_spectrogram(x_probe, fs, win_size_smp, win_shft_smp, [], [], name);
        h_tf = subplot(4,2,[3,4]);
        plot_spectrogram_data(x_spectrogram, [], [], freq_range, freq_color)
        set(gca, 'XLim', time_lim);

        % update plots
        update_plots();
    end

    % update plots when a new segment is selected
    function update_plots()
        % constrain ix_segment to [1 : n_segments]
        ix_segment = max([1, ix_segment]);
        ix_segment = min([n_segments, ix_segment]);

        % delete old areas (if existent)
        if exist('h_area1', 'var')
            delete(h_area1);
            delete(h_area2);
        end

        % select segment
        x = x_segments(:, :, ix_segment);

        % compute and plot Modulation Spectrogram
        clc;
        disp('computing modulation spectrogram...');
        x_stft_modspec = strfft_modulation_spectrogram(x, fs, win_size_smp, win_shft_smp, 2, [], 2, [], name);
        subplot(4,2,[6,8])
        plot_modulation_spectrogram_data(x_stft_modspec, [], freq_range, mfreq_range, mfreq_color)
        % Uncomment for log axes
        %set(gca,'XScale','log');
        %set(gca,'YScale','log');

        % plot spectrogram for segment
        subplot(4,2,7)
        plot_spectrogram_data(x_stft_modspec.spectrogram_data, [], [], freq_range, freq_color)

        % plot time series for segment
        subplot(4,2,5)
        plot_signal(x, fs, name);
        colorbar();

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
        fprintf('signal name            : %s\n', name{1} );
        fprintf('segment size  (seconds): %0.3f\n', seg_size_sec);
        fprintf('segment shift (seconds): %0.3f\n', seg_shft_sec);
        fprintf('segment position  (sec): %0.3f\n', seg_ini_sec)
        fprintf('window size   (seconds): %0.3f\n', win_size_sec);
        fprintf('window shift  (seconds): %0.3f\n', win_shft_sec);
        fprintf('windows per segment    : %d\n', x_stft_modspec.n_windows);

        drawnow();
    end

    % update parameters for AM analysis
    function update_parameters()
        % GUI to get new parameters
        prompt = {'Segment         (seconds): ', ...
                  'Segment shift   (seconds): ', ...
                  'Window size     (seconds): ', ...
                  'Window shift    (seconds): ', ...
                  'Freq Conv. min,Max (Hz): ', ...
                  'Spectr Pwr min,Max (dB): ', ...
                  'Freq Mod.  min,Max (Hz): ', ...
                  'ModSpec Pwr min,Max(dB): ', ...
                  };

        win_name = 'Modulation Analysis Parameters'; numlines = 1;

        defaultanswer = {num2str(seg_size_sec), num2str(seg_shft_sec), ...
                         num2str(win_size_sec), num2str(win_shft_sec), ...
                         num2str(freq_range), num2str(freq_color), ...
                         num2str(mfreq_range), num2str(mfreq_color),...
                         };

        options.Resize = 'off';

        answer = inputdlg(prompt,win_name,numlines,defaultanswer,options);

        % if the user clicks the 'Cancel' button, 'answer' is empty
        if ~isempty(answer)
            seg_size_sec = str2double(answer{1});  % (seconds)
            seg_shft_sec = str2double(answer{2});  % (seconds)
            win_size_sec = str2double(answer{3});  % (seconds)
            win_shft_sec = str2double(answer{4});  % (seconds)
            freq_range   = str2double(answer{5});  % (Hz)
            freq_color   = str2double(answer{6});  % (dB)
            mfreq_range  = str2double(answer{7});  % (Hz)
            mfreq_color  = str2double(answer{8});  % (dB)

            % call first_run()
            first_run();
        end
    end

% if the user press ESC, close the GUI and clean console
close(h_fig);
clc;

end
