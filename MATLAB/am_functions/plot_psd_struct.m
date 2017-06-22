function plot_psd_struct(psd_struct, ix, a_range, f_range)
% plot_psd_struct(spectrogram_struct, ix, a_range, f_range)
% Plots teh Spectrogram(s) in the spectrogram struct
%
% INPUTS:
%  psd_struct           = Spectrogram structure
%  ix                   = Array indicating the channels to be plotted
%                         (all channels if empty)
%                         If ix has one element, the plot will be
%                         located in the current axes, otherwise, a new
%                         figure will be created per plot
%  a_range              = Amplitude axis range [min, max] in (seconds)
%  f_range              = Frequency axis range [min, max] in (Hz)

% Validate 'type' argumet
if ~exist('ix','var') || isempty(ix)
    ix = 1 : size(psd_struct.PSD, 3);
end

% Check if ix has ONLY one element
if numel(ix) == 1
    new_figure = false;
    % Retrieve Current Axes handle from the Current Figure, if there is not
    % Current Figure, it's generated here
    h_ax = get(gcf, 'CurrentAxes');
    % If there is not Current Axes handle, create new Axes in current Figure
    if size(h_ax, 1) == 0
        h_ax = axes();
    end
else
    new_figure = true;
end

% Validate 'a_range' argumet
if ~exist('a_range','var')
    a_range = [];
end

% Validate 'f_range' argumet
if ~exist('f_range','var')
    f_range = [];
end

for i_channel = ix      
    if new_figure
            figure();
            h_ax = axes();
    end   
    plot_one_psd(h_ax, ...
                 psd_struct.PSD(:, :, i_channel), ...
                 psd_struct.freq_axis, ...
                 psd_struct.channel_names{i_channel},...
                 a_range, f_range);
end

    function plot_one_psd(h_ax, X_pwr, f_ax, title_str, a_range, f_range) 
        X_plot = 10*log10(X_pwr(:,:,1)' + eps);
        plot(h_ax, f_ax, X_plot );
        xlabel('Frequency (Hz)')
        ylabel('Power (dB/Hz)')
        set(h_ax, 'XMinorTick','on');
        set(h_ax, 'YMinorTick','on');
        if not(isempty(f_range))
            xlim(f_range);
        end
        if not(isempty(a_range))
            ylim(a_range);
        end
        title(title_str);
    end
end

