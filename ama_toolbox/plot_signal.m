function plot_signal(x, fs, name)
%PLOT_SIGNAL Behaves as matplotlib.pyplot.plot(x) but X axis is definded by `fs` [Hz]
%     
%     Parameters
%     ----------
%     x : 
%         1D Signal as column or row Vector 
%     fs :
%         Sampling frequency in Hz
%     name :
%         Name of the signal (Default 'Signal-01')

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

