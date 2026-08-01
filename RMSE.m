% Example data: Observed vs. Model Predicted values (e.g., log10 units/mL)
observed = [1.25, 2.40, 3.15, 2.80, 1.95];
predicted = [1.23, 2.42, 3.14, 2.82, 1.93];

% 1. Calculate Average Residual Difference (Mean Error)
mean_difference = mean(predicted - observed);

% 2. Calculate Root Mean Squared Error (RMSE)
% Method A: Using immse
rmse_a = sqrt(immse(predicted, observed));

% Method B: Manual calculation
rmse_b = sqrt(mean((predicted - observed) .^ 2));

fprintf('Average Difference: %.3f log10 units/mL\n', mean_difference);
fprintf('RMSE: %.3f log10 units/mL\n', rmse_b);