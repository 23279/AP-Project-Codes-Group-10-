%% 
%  Parameter Sensitivity antibody decay rate (r4)
clear; clc; close all;

% System Parameters & Initial Conditions
t_span = [0 45];
a0 = 0.1;
v0 = 0;
y0 = [a0; v0];

k1 = 10; k2 = 50; alpha = 1;
r1 = 0.019; r2 = 0.26; r3 = 0.014;
r4_base = 0.038; % Baseline decay rate

pulse = @(t) double(t>=0 & t<=1) + double(t>=21 & t<=22);

% Base ODE Function
get_ode = @(r4_val) @(t, y) [ ...
    r1*y(2) + r2*y(1)*y(2) + y(1)*(r3 - r4_val*y(1)); ...
    alpha*pulse(t) - (k1*y(2)) / (k2 + y(2) + eps) ...
];

% Run Simulations for the Three Scenarios

% Scenario A: Low Decay (Half)
[t_low, y_low] = ode45(get_ode(r4_base * 0.5), t_span, y0);

% Scenario B: Baseline
[t_base, y_base] = ode45(get_ode(r4_base), t_span, y0);

% Scenario C: High Decay (Double)
[t_high, y_high] = ode45(get_ode(r4_base * 2), t_span, y0);

% plot antibody trajectories
figure(9);
plot(t_low, y_low(:,1), 'g-', 'LineWidth', 1.5); hold on;
plot(t_base, y_base(:,1), 'b-', 'LineWidth', 1.5);
plot(t_high, y_high(:,1), 'r-', 'LineWidth', 1.5);

grid on;
title('parameter sensitivity: antibody decay rate (r4)');
xlabel('days');
ylabel('antibody level a(t)');
legend('low decay (r4 * 0.5)', 'baseline', 'high decay (r4 * 2)', 'Location', 'northeast');