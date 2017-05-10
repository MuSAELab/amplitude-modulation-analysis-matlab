function Y = rfft(X, n, dim)
%Y = RFFT(X, n, dim)
%
%RFFT Real Fast Fourier Transform
%   The RFFT function returns the DFT of a real signal (or signals) in X,
%   using the FFT algorithm.
%   The output Y contains ONLY the positive frequencies
%
%   The usage for RFFT() is the same as the one for FFT()
%   Refer to FFT() help for further details
%
%   Considering a real signal A with B = fft(A), B is Hermitian symmetric,
%   i.e. B(-1) = conj(B(1)), therefore the complete spectrum B
%   can be found by using with only the positive frequencies in B
%
%   The RFFT uses the Hermitian Symmetry property, returning ONLY positive frequencies
%   In this way the spectrum size has approx. 50% the size of the complete
%   spectrum, improving storage without losing information
%
% Example:
%
% n = size(X ,1); %Assuming one or more column-vector signals
% Y = RFFT(X);
% X_recover = IRFFT(Y, n);
%
% Raymundo Cassani

% verify X
shape_X = size(X);
if numel(shape_X) > 2
    error('RFFT only accepts 1D or 2D arrays as X input');
end

% check shape of X, and set n and dim defaults
if isvector(X)
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
Yc = fft(X, n, dim);

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
    Y = Yc(1:n_new, :);
else
    Y = Yc(:, 1: n_new);
end

end
