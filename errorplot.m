%% 
% Solver Error Analysis - Global Error Convergence
clear; clc; close all;

% System Parameters & Initial Conditions
t_span = [0 45];
a0 = 0.1; 
v0 = 0;   
y0 = [a0; v0]; 

k1 = 10; k2 = 50; alpha = 1;
r1 = 0.019; r2 = 0.26; r3 = 0.014; r4 = 0.038;

% Vaccination Pulse (Day 0 and Day 21)
pulse = @(t) double(t>=0 & t<=1) + double(t>=21 & t<=22);

% ODE Function
system_ode = @(t, y) [ ...
    r1*y(2) + r2*y(1)*y(2) + y(1)*(r3 - r4*y(1)); ...
    alpha*pulse(t) - (k1*y(2)) / (k2 + y(2) + eps) ...
];

%  Calculate Exact Reference Solution (ode45 with tight tolerances)
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11, 'MaxStep', 0.01);
[t_ref, y_ref] = ode45(system_ode, t_span, y0, options);

% Step Sizes 
dt_values = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01];
num_steps = length(dt_values);

err_euler = zeros(1, num_steps);
err_rk4 = zeros(1, num_steps);

% 4. Run Solvers and Calculate Error
for i = 1:num_steps
    dt = dt_values(i);
    t_grid = t_span(1):dt:t_span(2);
    N = length(t_grid);
    
    y_euler = zeros(2, N); y_euler(:,1) = y0;
    y_rk4 = zeros(2, N);   y_rk4(:,1) = y0;
    
    for n = 1:(N-1)
        t_curr = t_grid(n);
        
        % Euler
        f1 = system_ode(t_curr, y_euler(:,n));
        y_euler(:,n+1) = y_euler(:,n) + dt * f1;
        
        % RK4
        k_1 = system_ode(t_curr, y_rk4(:,n));
        k_2 = system_ode(t_curr + 0.5*dt, y_rk4(:,n) + 0.5*dt*k_1);
        k_3 = system_ode(t_curr + 0.5*dt, y_rk4(:,n) + 0.5*dt*k_2);
        k_4 = system_ode(t_curr + dt, y_rk4(:,n) + dt*k_3);
        y_rk4(:,n+1) = y_rk4(:,n) + (dt/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
    end
    
    % Interpolate Reference to find Max Global Error
    A_ref_interp = interp1(t_ref, y_ref(:,1), t_grid);
    err_euler(i) = max(abs(y_euler(1,:) - A_ref_interp));
    err_rk4(i)   = max(abs(y_rk4(1,:)   - A_ref_interp));
end


figure(1);
set(gcf, 'Color', 'w');

% Plot actual errors 
loglog(dt_values, err_euler, '-ob', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'none'); hold on;
loglog(dt_values, err_rk4, '-sr', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'none');


loglog(dt_values, 6 * dt_values.^1, 'k--', 'LineWidth', 1);      
loglog(dt_values, 1000 * dt_values.^4, 'k:', 'LineWidth', 1);     

set(gca, 'XDir', 'reverse', 'XScale', 'log', 'YScale', 'log');
xlim([1e-2, 1e0]);
ylim([1e-5, 1e2]);
grid on;
grid minor; 
title('Solver Error Analysis: Global Error Convergence', 'FontWeight', 'bold');
xlabel('Step Size dt (Days)');
ylabel('Max Absolute Global Error');


legend({'Euler Method (Actual)', 'RK4 Method (Actual)', ...
    'Theoretical \mathcal{O}(dt^1) Linear', 'Theoretical \mathcal{O}(dt^4) Quartic'}, ...
    'Location', 'southeast', 'Interpreter', 'tex');