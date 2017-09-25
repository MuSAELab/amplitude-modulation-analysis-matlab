function h = varea(X, color_cell, alpha_v)
%VAREA Summary of this function goes here
%   Detailed explanation goes here

if size(X, 2) ~= 2
    error('X must be of the size [n, 2]')
end

% Validate 'alpha_v' argumet
if ~exist('alpha_v','var') || isempty(alpha_v)
    alpha_v = 0.2;
end;

y_limits = get(gca,'ylim');
y_limits2 = [y_limits, fliplr(y_limits)];

n_areas = size(X, 1);

for i_area = 1 : n_areas
    x1 = X(i_area, 1); 
    x2 = X(i_area, 2);
    hold on;
    h = fill([x1, x1, x2, x2], y_limits2, color_cell{i_area});
    if ~is_octave()
        alpha(alpha_v);
    end
    hold off 
end


end

