

%  Variables

t = [0 45];  
a0 = 0.1;  
v0 = 0;     

% Added baseline parameters 
k1 = 10; 
k2 = 50; 
alpha = 1; 
r1 = 0.019; 
r2 = 0.26; 
r3 = 0.014; 
r4 = 0.038; 

% Number of distinct parameter profiles to simulate
num_samples = 500; 

% Explore a realistic +/- 20% variance around baseline values
range_pct = 0.20;  

param_names = {'k1', 'k2', 'alpha', 'r1', 'r2', 'r3', 'r4'};
base_p = [k1, k2, alpha, r1, r2, r3, r4];
num_p = length(base_p);

% Define parameter boundaries
lower_bounds = base_p * (1 - range_pct);
upper_bounds = base_p * (1 + range_pct);

% Generate Latin Hypercube Sample Matrix

rng(42); % Set random seed 
lhs_normalized = lhsdesign(num_samples, num_p);

% Scale normalized LHS matrix to our actual parameter bounds
LHS_matrix = lower_bounds + lhs_normalized .* (upper_bounds - lower_bounds);

% Run Monte Carlo Simulation Ensemble
simulated_peaks = zeros(num_samples, 1);

% Common pulse check expression inside loop
pulse_days = @(t_curr) ((t_curr >= 0.0 & t_curr <= 1.0) | (t_curr >= 21.0 & t_curr <= 22.0));

for s = 1:num_samples
    % Extract parameter set for this specific simulation sample
    sp = LHS_matrix(s, :);
    
    % Temporary ODE using sample parameter profile
    gsa_ode = @(t_curr, y) [ ...
        sp(4)*y(2) + sp(5)*y(1)*y(2) + y(1)*(sp(6) - sp(7)*y(1)); ... % dadt
        sp(3)*pulse_days(t_curr) - (sp(1)*y(2))/(sp(2) + y(2) + eps) ...   % dvdt added eps just in case to avoid dividing by zero
    ];
    
    [~, y_sim] = ode45(gsa_ode, t, [a0; v0]);
    simulated_peaks(s) = max(y_sim(:, 1)); % Save target output metric
end

% Compute Spearman Rank Correlation Coefficients

global_correlations = zeros(1, num_p);
for p = 1:num_p
    global_correlations(p) = corr(LHS_matrix(:, p), simulated_peaks, 'Type', 'Spearman');
end

%  Global Parameter Importance Rank

figure(13);
set(gcf, 'Position', [120, 120, 900, 500]);

h_gsa = bar(global_correlations, 'FaceColor', 'flat');

%  Green for positive drivers, Red for suppressors
for b = 1:num_p
    if global_correlations(b) >= 0
        h_gsa.CData(b,:) = [0.15, 0.55, 0.35]; 
    else
        h_gsa.CData(b,:) = [0.75, 0.25, 0.25]; 
    end
end

set(gca, 'XTick', 1:num_p, 'XTickLabel', param_names);
title('LHS Monte Carlo Sweep');
xlabel('Differential Equation Input Parameters');
ylabel('Spearman Rank Correlation Coefficient (\rho)');
grid on;
ylim([-1.1, 1.1]);

% Add value readouts directly above/below the bars
for b = 1:num_p
    val = global_correlations(b);
    if val >= 0
        text(b, val + 0.05, num2str(val, '%0.3f'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    else
        text(b, val - 0.08, num2str(val, '%0.3f'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end