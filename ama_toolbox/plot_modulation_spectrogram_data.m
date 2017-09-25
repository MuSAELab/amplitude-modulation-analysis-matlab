function plot_modulation_spectrogram_data(modulation_spectrogram_data, ix, f_range, fm_range, c_range, c_map)
% PLOT_MODULATION_SPECTROGRAM(modulation_spectrogram_data, ix, f_range, fm_range, c_range)
% Plot the Power Modulation Spectrogram related to the `modulation_spectrogram_data`
%         
%     Parameters
%     ----------
%     modulation_spectrogram_data : 
%         Structure with Modulation Spectrogram data
%     ix : Index of the signal (channel) to plot
%         (Default, all the channels, a new figure for each)
%     f_range : Frequency range
%         (Default [minimum frequency, maximum frequency])
%     fm_range : Modulation frequency range
%         (Default [minimum mod_frequency, maximum mod_frequency])
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
    ix = 1 : size (modulation_spectrogram_data.power_modulation_spectrogram, 3);
end

% Validate 'c_map' argumet
if ~exist('c_map','var') || isempty(ix)
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

% Validate 'f_range' argumet
if ~exist('f_range','var')
    f_range = [];
end

% Validate 'fm_range' argumet
if ~exist('fm_range','var')
    fm_range = [];
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
    plot_one_modspectrogram(h_ax, ...
                            modulation_spectrogram_data.power_modulation_spectrogram(:, :, i_channel), ...
                            modulation_spectrogram_data.freq_axis, ...
                            modulation_spectrogram_data.freq_mod_axis, ...
                            modulation_spectrogram_data.channel_names{i_channel}, ...
                            f_range, fm_range, c_range, c_map);
end

    function plot_one_modspectrogram(h_ax, X_pwr, f_ax, fmod_ax, title_str, f_range, fm_range, c_range, c_map)      
        X_plot = double(10*log10(X_pwr(:,:,1) + eps));
        surf(h_ax, fmod_ax, f_ax, X_plot,'LineStyle','none')
        colormap(c_map);
        view(0,90)
        xlabel('modulation frequency (Hz)');
        ylabel('conventional frequency (Hz)');
        set(h_ax, 'XMinorTick','on');
        set(h_ax, 'YMinorTick','on');
        if not(isempty(fm_range))
            xlim(fm_range);
        else
            xlim([min(fmod_ax), max(fmod_ax)]);
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

