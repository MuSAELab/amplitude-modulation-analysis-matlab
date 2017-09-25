function plot_psd_data(psd_data, ix, a_range, f_range)
% PLOT_PSD_DATA(spectrogram_struct, ix, a_range, f_range)
% Plot the PSD related to the `psd_data`
%         
%     Parameters
%     ----------
%     psd_data : 
%         Structure with PSD data
%     ix : Index of the signal (channel) to plot
%         (Default, all the channels, a new figure for each)
%     p_range : Power range
%         (Default [minimum power, maximum power])
%     f_range : Frequency range
%         (Default [minimum frequency, maximum frequency])
%     
%     Returns
%     -------
%     If only a plot is requested, it is plotted in the existen axes (created if needed)
%     If many plots are requested, a new figure is created for each plot

% Validate 'type' argumet
if ~exist('ix','var') || isempty(ix)
    ix = 1 : size(psd_data.PSD, 3);
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
                 psd_data.PSD(:, :, i_channel), ...
                 psd_data.freq_axis, ...
                 psd_data.channel_names{i_channel},...
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

