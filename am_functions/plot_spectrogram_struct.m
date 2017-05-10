function plot_spectrogram_struct(spectrogram_struct, ix, t_range, f_range, c_range)
% plot_spectrogram_struct(spectrogram_struct, type)
% Plots teh Spectrogram(s) in the spectrogram struct
%
% INPUTS:
%  spectrogram_struct   = Spectrogram structure
%  ix                   = Array indicating the channels to be plotted
%                         (all channels if empty)
%                         If ix has one element, the plot will be
%                         located in the current axes, otherwise, a new
%                         figure will be created per plot
%  t_range              = Time axis range [min, max] in (seconds)
%  f_range              = Frequency axis range [min, max] in (Hz)
%  c_range              = Color axis range [min, max] in (dB)     

% Validate 'type' argumet
if ~exist('ix','var') || isempty(ix)
    ix = 1 : size(spectrogram_struct.pwr_spectrogram, 3);
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
                         spectrogram_struct.pwr_spectrogram(:, :, i_channel), ...
                         spectrogram_struct.time_axis, ...
                         spectrogram_struct.freq_axis, ...
                         spectrogram_struct.channel_names{i_channel},...
                         t_range, f_range, c_range);
end

    function plot_one_spectrogram(h_ax, X_pwr, t_ax, f_ax, title_str, t_range, f_range, c_range) 
        X_plot = double(10*log10(X_pwr(:,:,1)' + eps));
        surf(h_ax, t_ax, f_ax, X_plot,'LineStyle','none')
        view(0,90)
        xlabel('time (s)');
        ylabel('frequency (Hz)');
        set(h_ax, 'XMinorTick','on');
        set(h_ax, 'YMinorTick','on');
        if not(isempty(t_range))
            xlim(t_range);
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

