function plot_modspectrogram_struct(modspectrogram_struct, ix, f_range, fm_range, c_range)
% plot_modspectrogram_struct(spectrogram_struct, ix, f_range, fm_range, c_range)
% Plots the Modulation Spectrogram(s) in the modspectrogram_struct
%
% INPUTS:
%  modspectrogram_struct   = Modulation Spectrogram structure
%  Optional:
%   ix                     = Array indicating the channels to be plotted
%                             (all channels if empty)
%                             If ix has one element, the plot will be
%                             located in the current axes, otherwise, a new
%                             figure will be created per plot
%   f_range                = Conventional frequency axis range [min, max] in (Hz)
%   fm_range               = Modulation frequency axis range [min, max] in (Hz)
%   c_range                = Color axis range [min, max] in (dB) 

% Validate 'ix' argumet
if ~exist('ix','var') || isempty(ix)
    ix = 1 : size (modspectrogram_struct.pwr_modspec, 3);
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
                            modspectrogram_struct.pwr_modspec(:, :, i_channel), ...
                            modspectrogram_struct.freq_axis, ...
                            modspectrogram_struct.modfreq_axis, ...
                            modspectrogram_struct.channel_names{i_channel}, ...
                            f_range, fm_range, c_range);
end

    function plot_one_modspectrogram(h_ax, X_pwr, f_ax, fmod_ax, title_str, f_range, fm_range, c_range)      
        X_plot = double(10*log10(X_pwr(:,:,1) + eps));
        surf(h_ax, fmod_ax, f_ax, X_plot,'LineStyle','none')
        view(0,90)
        xlabel('modulation frequency (Hz)');
        ylabel('conventional frequency (Hz)');
        set(h_ax, 'XMinorTick','on');
        set(h_ax, 'YMinorTick','on');
        if not(isempty(fm_range))
            xlim(fm_range);
        end
        if not(isempty(f_range))
            ylim(f_range);
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

