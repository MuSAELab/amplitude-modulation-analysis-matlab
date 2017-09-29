function plot_spectrogram_data(spectrogram_data, ix, t_range, f_range, c_range, c_map)
% PLOT_SPECTROGRAM_DATA(spectrogram_data, ix, t_range, f_range, c_range, c_map)
% Plot the Power Spectrogram related to the `spectrogram_data`
%         
%     Parameters
%     ----------
%     spectrogram_data : 
%         Structure with Spectrogram data
%     ix : Index of the signal (channel) to plot
%         (Default, all the channels, a new figure for each)
%     t_range : Time range
%         (Default [minimum time, maximum time])
%     f_range : Frequency range
%         (Default [minimum frequency, maximum frequency])
%     c_range : Color (power) range
%         (Default [mean power, maximum power])
%     c_map : Colot Map
%         (Default viridis)
%    
%     Returns
%     -------
%     If only a plot is requested, it is plotted in the existen axes (created if needed)
%     If many plots are requested, a new figure is created for each plot     

% Validate 'ix' argumet
if ~exist('ix','var') || isempty(ix)
    ix = 1 : size(spectrogram_data.power_spectrogram, 3);
end

% Validate 'c_map' argumet
if ~exist('c_map','var') || isempty(c_map)
    c_map = 'viridis';
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

% Validate 't_range' argumet
if ~exist('t_range','var')
    t_range = [];
end

% Validate 'f_range' argumet
if ~exist('f_range','var')
    f_range = [];
end

% Validate 'c_range' argumet
if ~exist('c_range','var')
    c_range = [];
end


for i_channel = ix      
    if new_figure
            figure();
            h_ax = axes();
    end   
    plot_one_spectrogram(h_ax, ...
                         spectrogram_data.power_spectrogram(:, :, i_channel), ...
                         spectrogram_data.time_axis, ...
                         spectrogram_data.freq_axis, ...
                         spectrogram_data.channel_names{i_channel},...
                         t_range, f_range, c_range, c_map);
end

    function plot_one_spectrogram(h_ax, X_pwr, t_ax, f_ax, title_str, t_range, f_range, c_range, c_map) 
        X_plot = double(10*log10(X_pwr(:,:,1)' + eps));
        surf(h_ax, t_ax, f_ax, X_plot,'LineStyle','none')
        colormap(c_map);
        view(0,90)
        xlabel('time (s)');
        ylabel('frequency (Hz)');
        set(h_ax, 'XMinorTick','on');
        set(h_ax, 'YMinorTick','on');
        if not(isempty(t_range))
            xlim(t_range);
        else
            xlim([min(t_ax), max(t_ax)]);
        end
        if not(isempty(f_range))
            ylim(f_range);
        else
            ylim([min(f_ax), max(f_ax)]);
        end
        if not(isempty(c_range))
            caxis(c_range);
        else
            caxis([mean(mean(X_plot)) max(max(X_plot))]);  
        end
        colorbar;       
        title(title_str);
    end
end

