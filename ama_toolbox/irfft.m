function X = irfft(Y, n, dim)
%x = IRFFT(y, n, dim) Inverse Real Fast Fourier Transform.   
%     The IRFFT function returns the Inverse DFT (using the RFFT algorithm)of
%     a spectrum Y containing ONLY the positive frequencies, with the
%     assumption than Y is the positive half of a Hermitian Symmetric spectrum
%     from a real signal X.
%         
%     Parameters
%     ----------
%     y : 1D or 2D array with the positive spectrum of 
%         real-valued signals with shape (n_samples, n_channels)
%     n : Number of samples in the original x signals 
%         N not provided. Y is assumed be obtained from a signal X with even number fo samples 
%     dim : Dimension to compute the IRFFT (Default: Last dimension in `y`) 
% 
%     Returns
%     -------
%     x : Real-valued signal(s) 
%     
%     
% See also IFFT
%
% Example:
%
% n = size(X ,1); %Assuming one or more column-vector signals
% Y = RFFT(X);
% X_recover = IRFFT(Y, n);
%


% verify Y
shape_Y = size(Y);
if numel(shape_Y) > 2
    error('IRFFT only accepts 1D or 2D arrays as Y input');
end

% check shape of Y, and set n and dim defaults
if isvector(Y)
    if shape_Y(1) == 1
        % Y is a row vector
        dim_def = 2;
    else
        % Y is a column vector
        dim_def = 1;
    end
else
    % Y is a 2D Matrix, a shape [n_samples, n_channels] is asummed
    dim_def = 1;
end

% verify 'dim' dimension parameter
if ~exist('dim','var') || isempty(dim)
    dim = dim_def;
end

% verify 'n' number-of-samples parameter
if ~exist('n','var') || isempty(n)
    warning('N not provided. Y is assumed be obtained from a signal X with even number fo samples');
    n_half = size(Y, dim);
    n = (n_half - 1) * 2;
end

% reconstruct missing half of Spectrum
if ~mod(n, 2)
    % number of samples is even
    n_half = (n / 2) + 1;
    ix_limit = (2 : n_half - 1 );
else
    % number of samples is odd
    n_half = (n + 1) / 2;
    ix_limit = (2 : n_half );
end

% check shape of Y, and add negative frequencies
if dim == 1
    % spectra in Y are column wise
    Y_neg = conj(flipud(Y(ix_limit, :)));
    Yc = [Y; Y_neg];
else
    % spectra in Y are row-wise
    Y_neg = conj(fliplr(Y(:, ix_limit)));
    Yc = [Y, Y_neg];
end

X = real(ifft(Yc, n, dim));

end
