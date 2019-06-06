function y = rfft(x, n, dim)
%y = RFFT(x, n, dim) Real Fast Fourier Transform.   
%     Considering a real signal A with B = fft(A), B is Hermitian symmetric,
%     i.e. B(-1) = conj(B(1)), therefore the complete spectrum B
%     can be found by using with only the non-negative frequencies in B
%         
%     Parameters
%     ----------
%     x : 1D array with shape (n_samples) or 2D array with shape (n_samples, n_channels)
%     n : Number of samples to compute the FFT
%     dim : Dimension to compute the RFFT (Default: Last dimension in `x`) 
% 
%     Returns
%     -------
%     y : Non-negative complex spectrum of `x`, with shape as `x` 
%     
%     
% See also FFT
%
% Example:
% xi = rand(100,1);
% xi_rfft = rfft(xi);
% xo = irfft(xi_rfft, 100);
%

% verify X
shape_X = size(x);
if numel(shape_X) > 2
    error('RFFT only accepts 1D or 2D arrays as X input');
end

% check shape of X, and set n and dim defaults
if isvector(x)
    if shape_X(1) == 1
        % X is a row vector
        dim_def = 2;
    else
        % X is a column vector
        dim_def = 1;
    end
else
    % X is a 2D Matrix, a shape [n_samples, n_channels] is asummed
    dim_def = 1;
end

% verify 'dim' dimension argument
if ~exist('dim','var') || isempty(dim)
    dim = dim_def;
end

% verify 'n'
if ~exist('n','var') || isempty(n)
    n = shape_X(dim);
end

% FFT
Yc = fft(x, n, dim);

% points to keep
if ~mod(n,2)
    % number of samples is even
    n_new = (n / 2) + 1;
else
    % number of samples is odd
    n_new = (n + 1) / 2;
end

% remove negative frequencies
if dim == 1
    y = Yc(1:n_new, :);
else
    y = Yc(:, 1: n_new);
end

end
