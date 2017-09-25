function c = conv_fft(a, b, mode)
%CONV_FFT Convolve the a vector with collection of vectors. 
%
%   c = CONV_FFT(a, b) 
%
%   Convolve a 1D array `a` with each column of the 2D array `b`.
%   `c` has the same class as `a` and `b`
%   This function is based in the Convolution Theorem.
%   c = ifft(fft(a).*fft(b(:,k)) for k = 1, ... K kernels
%
%   Paramters
%       a    :  1D array
%       b    :  2D array
%       mode :  string 'full' (default) | 'same'
%               'full'  - The output is the full discrete linear convolution
%                         of the inputs. (Default)
%               'same'  - The output is the same size as `a`, centered
%                         with respect to the 'full' output.
%   Returns
%       c : A 2D array where each columns corresponds to the 
%           convolution of `a` and a column of `b`
%
%   See also DECONV, CONV2, CONV, FILTER and,
%   in the Signal Processing Toolbox, XCORR, CONVMTX.


if ~isvector(a) || ~ismatrix(b)
  error(message('MATLAB:conv:ANotVector or BNotMatrix'));
end

% Validate 'shape' argument
if ~exist('mode','var') || isempty(mode)
    mode = 'full';
end

if ~ischar(mode)
  error(message('MATLAB:conv:unknownShapeParameter'));
end

% Input A vector to row vector
a = a(:);
sza = numel(a);

% Full convolution is computed
N = sza + size(b,1) -1;

ffta = fft(a,N);

K = size(b,2);
ctmp = zeros(N,K);

if K == 1
    ctmp = ifft(ffta.*fft(b,N),N);
else
    for ik=1:K
        ctmp(:,ik) = ifft(ffta.*fft(b(:,ik),N),N);
    end
end

if mode(1) == 'f'
   c = ctmp;
    return
end

if mode(1) == 's'
    tocut = (N-sza)/2;
    if K==1
        c = ctmp(1+ceil(tocut):end-floor(tocut));
        return
    else
        c = ctmp(1+ceil(tocut):end-floor(tocut),:);
        return
    end
end

end
