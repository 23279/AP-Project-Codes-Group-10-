% --- Graph 1: Residual Plot ---
observed = [1.25, 2.40, 3.15, 2.80, 1.95];
predicted = [1.23, 2.42, 3.14, 2.82, 1.93];
time_points = [0, 7, 14, 21, 28]; % Example time points (days)

residuals = predicted - observed;

figure;
stem(time_points, residuals, 'filled', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
yline(0, '--r', 'LineWidth', 1.2); % Zero-error reference line
hold off;

xlabel('Time (Days)');
ylabel('Residuals (Predicted - Observed)');
title('Model Residual Plot');
grid on;