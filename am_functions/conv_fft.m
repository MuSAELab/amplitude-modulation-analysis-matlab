function c = conv_fft(a, b, shape)
%CONV Convolution and polynomial multiplication.
%   C = CONV_FFT(A, B) convolves vectors A and a vector or matrix B;
%   if B is a vector, C is the convolution between A and B
%   if B is a matrix (B is a collection of K kernels), [samples,K]
%   C is matrix that  contains the result of the convolution of A
%   and each kernel in B. C = [convolution_length, K]
%
%  C has the same class as A and B
%
%   This function is based in the Convolution Theorem.
%   C = ifft(fft(A).*fft(B(:,k)) for k = 1, ... K kernels
%
%   The resulting convolution vector is
%   length MAX([LENGTH(A)+LENGTH(B)-1,LENGTH(A),LENGTH(B)]).
%   If A and B are vectors of polynomial coefficients, convolving
%   them is equivalent to  multiplying the two polynomials.
%
%   C = CONV_FFT(A, B, SHAPE) returns a subsection of the convolution with size
%   specified by SHAPE:
%     'full'  - (default) returns the full convolution,
%     'same'  - returns the central part of the convolution
%               that is the same size as A.
%
%   See also DECONV, CONV2, CONV, FILTER and,
%   in the Signal Processing Toolbox, XCORR, CONVMTX.
%
%   Raymundo Cassani.
%   June 2014

if ~isvector(a) || ~ismatrix(b)
  error(message('MATLAB:conv:ANotVector or BNotMatrix'));
end

% Validate 'shape' argument
if ~exist('shape','var') || isempty(shape)
    shape = 'full';
end

if ~ischar(shape)
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

if shape(1) == 'f'
   c = ctmp;
    return
end

if shape(1) == 's'
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
