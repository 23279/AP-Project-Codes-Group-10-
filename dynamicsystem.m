%% 
% Dynamic system structural eigenvalues & Instability
clear; clc; close all;

% 1. System Parameters & Initial Conditions
t_span = [0 45];
y0 = [0.1; 0]; % [A(0); V(0)]

k1 = 10; k2 = 50; alpha = 1;
r1 = 0.019; r2 = 0.26; r3 = 0.014; r4 = 0.038;

% Vaccination Pulse (Day 0 and Day 21)
pulse = @(t) double(t>=0 & t<=1) + double(t>=21 & t<=22);

% ODE System
sys_ode = @(t,y) [ ...
    r1*y(2) + r2*y(1)*y(2) + y(1)*(r3 - r4*y(1)); ...
    alpha*pulse(t) - (k1*y(2)) / (k2 + y(2) + eps) ...
];

% --- Subplot 1: Structural Eigenvalues ---
[t_ref, y_ref] = ode45(sys_ode, t_span, y0);
A = y_ref(:,1);
V = y_ref(:,2);

% Calculate Jacobian diagonal elements over time
lambda_1 = r2.*V + r3 - 2.*r4.*A;          % Antibody Dynamics
lambda_2 = -(k1.*k2) ./ ((k2 + V).^2);     % Cell Decay Speed

figure(7);
set(gcf, 'Position', [100, 100, 1000, 450], 'Color', 'w');

subplot(1,2,1);
plot(t_ref, lambda_1, 'LineWidth', 1.5); hold on;
plot(t_ref, lambda_2, 'LineWidth', 1.5);
title('Dynamic System Structural Eigenvalues', 'FontWeight', 'bold');
xlabel('Days'); 
ylabel('Eigenvalue Real Magnitude (\lambda)');
legend('\lambda_1 (Antibody Dynamics)', '\lambda_2 (Cell Decay Speed)', 'Location', 'best');
grid on; grid minor;
ylim([-0.2, 0.4]);

% Euler Instability at dt = 6.793 
dt_unstable = 6.793;
t_euler = t_span(1):dt_unstable:t_span(2);
num_steps = length(t_euler);
y_euler = zeros(2, num_steps);
y_euler(:,1) = y0;

for n = 1:(num_steps-1)
    f = sys_ode(t_euler(n), y_euler(:,n));
    y_euler(:,n+1) = y_euler(:,n) + dt_unstable * f;
end

subplot(1,2,2);
% Plot the stable ode45 reference and unstable Euler 
plot(t_ref, y_ref(:,1), 'Color', [0.5 0.6 0.2], 'LineWidth', 2); hold on;
plot(t_euler, y_euler(1,:), '-or', 'LineWidth', 1.2, 'MarkerFaceColor', 'r');

title('Numerical Instability at dt = 6.793 Days', 'FontWeight', 'bold');
xlabel('Days'); 
ylabel('Antibody Units a(t)');
legend('ode45 Reference (Stable)', 'Euler Method (Unstable)', 'Location', 'northwest');
grid on; grid minor;
xlim([0, 45]);
ylim([0, 0.55]);