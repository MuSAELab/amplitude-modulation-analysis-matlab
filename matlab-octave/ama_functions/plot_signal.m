function plot_signal(x, fs, name)
%PLOT_SIGNAL Behaves as plot(x) but X axis is definded by FS [Hz]
%
% X     1D Signal as column or row Vector 
% FS    Sampling frequency in Hz
% NAME  Name of the signal (Default 'Sigbnal-01')

% Create time vector
time_vector = (0:size(x,1)-1)./fs;

plot(time_vector,x)
xlabel('Time (s)')

if ~exist('name', 'var') || isempty(name);
    name = 'Signal-01';
end

set(gca, 'XLim', [0, time_vector(end)]);

title(name);

end

