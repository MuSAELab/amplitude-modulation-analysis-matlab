function x = iepoching(epochs, shift_epoch)
% x = IEPOCHING(epochs, shift_epochs)
%
%     Merges a set of epochs [n_samples_epoch, n_channels] into  
%     the complete signal(s) x [n_samples, n_channels] taking into account
%     the shift between consecutive epochs
%    
%     Parameters
%     ----------
%     epochs : 2D array_like with shape (n_samples, n_channels)
%     shift_epoch : number of samples in smaller epochs
% 
%     Returns
%     -------
%     x : 2D array with shape (samples_epoch, n_channels, n_epochs)

% Verify the number of Input arguments
if nargin < 2
    error('iepoching funtion requires 2 arguments');
end

% Obtain parameters
[size_epoch, n_channels, n_epochs] = size(epochs);
n_samples = (shift_epoch * (n_epochs - 1)) + size_epoch;
ix = 1 + [0:n_epochs - 1 ] * shift_epoch;

% merging matrix
merging = zeros(n_samples, n_channels, 2);
% Number of epochs that contribute for a specific point
n_merging = merging;

for i_epoch = 1 : n_epochs
    merging(ix(i_epoch) : ix(i_epoch) + size_epoch - 1, :, 2 ) = squeeze(epochs(:, :, i_epoch)); 
    n_merging(ix(i_epoch) : ix(i_epoch) + size_epoch - 1, :, 2) = 1;
    merging(:,:,1) = sum(merging, 3);
    n_merging(:,:,1) = sum(n_merging, 3);
    merging(ix(i_epoch) : ix(i_epoch) + size_epoch - 1, :, 2 ) = 0;
    n_merging(ix(i_epoch) : ix(i_epoch) + size_epoch - 1, :, 2 ) = 0;
end

x = merging(:,:,1)./n_merging(:,:,1);