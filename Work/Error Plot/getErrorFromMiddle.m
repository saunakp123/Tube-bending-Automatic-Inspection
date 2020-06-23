function [lens, angle_errors, fit_errors] = getErrorFromMiddle(edge, true_angle)
    n = length(edge);
    lens = [];
    angles = [];
    fit_errors = [];
    % Compute the mid point
    mid = ceil((n + 1) / 2);
    for i = mid+5:5:n
        section = edge(2*mid-i:i, :);
        dis = norm(section(1, :) - section(end, :));
        p = polyfit(section(:, 2), section(:, 1), 1);
        angle = atand(p(1));
        angle = 90 + angle;
        angles = [angles, angle];

        f = polyval(p,section(:, 2));
        lens = [lens, dis];
        % Compute the fit error
        error = f - section(:, 1);
        error = floor(error);
        error = sqrt(sum(error.^2)) / length(error);
        fit_errors = [fit_errors, error];
    end
    angle_errors = abs(angles - true_angle);
end
