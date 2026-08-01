% global parameters
k1 = 10.0;
k2 = 50.0;
alpha = 1.0;

% initial conditions at day 0
a0 = 0.1;
v0 = 0.0;

% pfizer data
% subject id, r1, r2, r3, r4
pfizer_data = [
    3, 0.0062, 0.40, 0.0072, 0.027;
    4, 0.0110, 0.35, 0.0079, 0.026;
    19, 0.1100, 0.04, 0.0064, 0.028;
    41, 0.0130, 0.26, 0.0063, 0.028;
    44, 0.0210, 0.18, 0.0061, 0.028;
    76, 0.0250, 0.21, 0.0061, 0.029;
    146, 0.0210, 0.26, 0.0069, 0.027;
    179, 0.0110, 0.32, 0.0110, 0.067;
    215, 0.0100, 0.21, 0.0053, 0.029;
    216, 0.0160, 0.19, 0.0053, 0.030;
    226, 0.0061, 0.49, 0.0069, 0.027
];

% time setup
t_start = 0;
t_end = 45;
dt = 0.1;
t = t_start:dt:t_end;
num_steps = length(t);

% pick parameters for the first pfizer patient for the test
r1 = pfizer_data(1, 2);
r2 = pfizer_data(1, 3);
r3 = pfizer_data(1, 4);
r4 = pfizer_data(1, 5);

% storage arrays for euler
a_euler = zeros(1, num_steps);
v_euler = zeros(1, num_steps);
a_euler(1) = a0;
v_euler(1) = v0;

% euler method loop with standard pulse check matching document text
for i = 1:(num_steps - 1)
    current_t = t(i);
    current_a = a_euler(i);
    current_v = v_euler(i);
    
    % discrete pulse function for vaccine days
    u_t = 0.0;
    if (current_t >= 0.0 && current_t <= 1.0) || (current_t >= 21.0 && current_t <= 22.0)
        u_t = 1.0;
    end
    
    % system differential equations
    dadt = r1 * current_v + r2 * current_a * current_v + current_a * (r3 - r4 * current_a);
    dvdt = alpha * u_t - (k1 * current_v) / (k2 + current_v);
    
    % calculate next values using euler step
    a_euler(i + 1) = current_a + dadt * dt;
    v_euler(i + 1) = current_v + dvdt * dt;
end

% runge-kutta loop setup
% helper function for dynamic pulse evaluation
evaluate_u = @(tau) double((tau >= 0.0 && tau <= 1.0) || (tau >= 21.0 && tau <= 22.0));

% storage arrays for rk4
a_rk4 = zeros(1, num_steps);
v_rk4 = zeros(1, num_steps);
a_rk4(1) = a0;
v_rk4(1) = v0;

for i = 1:(num_steps - 1)
    t_now = t(i);
    a_now = a_rk4(i);
    v_now = v_rk4(i);
    
    % stage 1 (euler slopes)
    u1 = evaluate_u(t_now);
    dadt1 = r1 * v_now + r2 * a_now * v_now + a_now * (r3 - r4 * a_now);
    dvdt1 = alpha * u1 - (k1 * v_now) / (k2 + v_now);
    
    % stage 2 (midpoint using stage 1)
    t_mid = t_now + dt/2;
    a_mid1 = a_now + dadt1 * (dt / 2);
    v_mid1 = v_now + dvdt1 * (dt / 2);
    u2 = evaluate_u(t_mid); 
    
    dadt2 = r1 * v_mid1 + r2 * a_mid1 * v_mid1 + a_mid1 * (r3 - r4 * a_mid1);
    dvdt2 = alpha * u2 - (k1 * v_mid1) / (k2 + v_mid1);
    
    % stage 3 (midpoint using stage 2)
    a_mid2 = a_now + dadt2 * (dt / 2);
    v_mid2 = v_now + dvdt2 * (dt / 2);
    u3 = evaluate_u(t_mid); 
    
    dadt3 = r1 * v_mid2 + r2 * a_mid2 * v_mid2 + a_mid2 * (r3 - r4 * a_mid2);
    dvdt3 = alpha * u3 - (k1 * v_mid2) / (k2 + v_mid2);
    
    % stage 4 (end of step using stage 3)
    t_next = t_now + dt;
    a_end = a_now + dadt3 * dt;
    v_end = v_now + dvdt3 * dt;
    u4 = evaluate_u(t_next); 
    
    dadt4 = r1 * v_end + r2 * a_end * v_end + a_end * (r3 - r4 * a_end);
    dvdt4 = alpha * u4 - (k1 * v_end) / (k2 + v_end);
    
    % combine all stages
    a_rk4(i + 1) = a_now + (dt / 6) * (dadt1 + 2*dadt2 + 2*dadt3 + dadt4);
    v_rk4(i + 1) = v_now + (dt / 6) * (dvdt1 + 2*dvdt2 + 2*dvdt3 + dvdt4);
end

% ode45 implementation
ode_system = @(t_curr, y) [ ...
    r1 * y(2) + r2 * y(1) * y(2) + y(1) * (r3 - r4 * y(1)); ...
    alpha * ((t_curr >= 0.0 & t_curr <= 1.0) | (t_curr >= 21.0 & t_curr <= 22.0)) - (k1 * y(2)) / (k2 + y(2)) ... 
];

% solve using ode45
[~, y_ode] = ode45(ode_system, t, [a0; v0]);
a_ode45 = y_ode(:, 1)';

% plotting comparison between the methods
% euler vs rk4
figure('Color', 'w', 'Name', 'Solver Comparison: Euler vs RK4');
plot(t, a_euler, 'r--', 'LineWidth', 2); 
hold on;
plot(t, a_rk4, 'b-', 'LineWidth', 1.5);
xlabel('Days'); ylabel('Antibody Level A(t)');
title('Numerical Methods: Euler vs. RK4');
legend('Euler Method', 'RK4 Method', 'Location', 'best');
grid on; 
hold off;

% euler vs ode45
figure('Color', 'w', 'Name', 'Solver Comparison: Euler vs ode45');
plot(t, a_euler, 'r--', 'LineWidth', 2); 
hold on;
plot(t, a_ode45, 'k-', 'LineWidth', 1.5);
xlabel('Days'); ylabel('Antibody Level A(t)');
title('Numerical Methods: Euler vs. ode45');
legend('Euler Method', 'ode45 Built-in', 'Location', 'best');
grid on; 
hold off;

% rk4 vs ode45
figure('Color', 'w', 'Name', 'Solver Comparison: RK4 vs ode45');
plot(t, a_rk4, 'b--', 'LineWidth', 2); 
hold on;
plot(t, a_ode45, 'k-', 'LineWidth', 1.5);
xlabel('Days'); ylabel('Antibody Level A(t)');
title('Numerical Methods: Custom RK4 vs. ode45');
legend('Custom RK4', 'ode45 Built-in', 'Location', 'best');
grid on; 
hold off;