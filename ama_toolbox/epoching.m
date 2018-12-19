function [epochs, remainder, ix_center] = epoching(data, samples_epoch, samples_overlap)
% [epochs, remainder, ix_center] = EPOCHING(data, size_epoch, overlap_epoch)
%     Divides the `data` provided as [n_samples, n_channels] using the 
%     `size_epoch` indicated (in samples) and the `overlap_epoch` between 
%     consecutive epochs.
%    
%     Parameters
%     ----------
%     data : 2D array_like with shape (n_samples, n_channels)
%     samples_epochs : number of samples in smaller epochs
%     samples_overlap : number of samples for ovelap between epochs (Default 0)
% 
% 
%     Returns
%     -------
%     epochs : 3D array with shape (samples_epoch, n_channels, n_epochs)
%     remainder : 2D array with the remaining data after last complete epoch
%     ix_center : 1D array, indicates the index tha corresponds to the center of the
%                 nth epoch.
% 
%
%
% e.g 
% samples_epoch = 5, samples_overlap = 0
%
%  n_epochs = floor(n_samples / samples_epoch) = floor (23 / 5) = 4
%   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3       (Data array)
%  |E1       |E2       |E3       |E4       |Remainder
%       *         *         *         *                 (* = center of the epoch)      
% ix_center =  [3, 8, 13, 18]
%
% e.g 
% samples_epoch = 5  and samples_overlap = 2
% samples_shift = samples_epoch - samples_overlap 
%
% n_epochs = floor( (n_samples - samples_epoch) / samples_shift ) + 1 ;
% n_epochs = floor( (24 - 5) / 3 ) + 1  = 7
% 
%  1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4       (Data array)
% |E1-------| |E3-------| |E5-------| |E7-------|
%       |E2-------| |E4-------| |E6-------|     |R| 
%      *           *           *           *           
%            *           *           *                   (* = center of the epoch) 
% ix_center = [3, 6, 9, 12, 15, 18, 21]
% remainder = [4]

% Get class of data
data_class = class(data);

% Verify the number of Input arguments
if nargin < 3
    samples_overlap = 0;
end

% Obtain parameters of the data
n_samples = size(data,1);
n_channels = size(data,2);

% Size of half epoch
half_epoch = ceil(samples_epoch / 2);

% Epoch shift
samples_shift = samples_epoch - samples_overlap;

% Number of epochs
n_epochs = floor( (n_samples - samples_epoch) / samples_shift ) + 1 ;

%markers indicates where the epoch starts, and the epoch contains 
%size_epoch elements

markers = ((1:n_epochs-1)'*samples_shift)+1;
markers = [1; markers];

%Divide data in epochs
epochs = zeros(samples_epoch,n_channels,n_epochs, data_class);
ix_center = zeros(n_epochs,1, data_class);

if n_epochs == 0
    remainder = data;
    return
end

for i_epoch = 1:n_epochs
    epochs(:,:,i_epoch) = data( markers(i_epoch) : markers(i_epoch) + samples_epoch -1 ,:);
    ix_center(i_epoch) = markers(i_epoch) -1 + half_epoch;
end

if (markers(end) + samples_epoch - 1 < n_samples) 
    remainder = data(markers(end) + samples_epoch : n_samples,:);
else
    remainder = [];
end