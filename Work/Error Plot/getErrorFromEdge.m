function [lens, angle_errors, fit_errors] = getErrorFromEdge(edge, true_angle)
    edge = edge(50:end, :);
    n = length(edge);
    lens = [];
    angles = [];
    fit_errors = [];
    for i = 10:10:n
        section = edge(1:i, :);
        dis = norm(section(1, :) - section(end, :));
        p = polyfit(section(:, 2), section(:, 1), 1);
        angle = atand(p(1));
        angle = 90 + angle;
        angles = [angles, angle];

        f = polyval(p,section(:, 2));
        lens = [lens, dis];
        error = f - section(:, 1);
        error = floor(error);
        error = sqrt(sum(error.^2)) / length(error);
        fit_errors = [fit_errors, error];
    end
    angle_errors = abs(angles - true_angle);
end