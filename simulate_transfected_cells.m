
function simulate_transfected_cells()
% initialize the variables
    k1 = 10.0; 
    k2 = 50.0; 
    alpha = 1.0; 
    
% Set time span of 45 days with a step size of 0.1
    tspan = 0:0.1:45; 
    
% Set initial conditions for A and V
    y0 = [0.1; 0]; 

% Doses at day 0 and 21 for Pfizer
    pfizer_days = [0.0, 21.0]; 
    
% (Subject_ID, r1, r2, r3, r4)
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

% Doses at day 0 and 28 for Moderna
    moderna_days = [0.0, 28.0]; 
    
    moderna_data = [
        1, 0.0190, 0.26, 0.0140, 0.038;
        5, 0.0450, 0.17, 0.0140, 0.025;
        20, 0.0260, 0.23, 0.0140, 0.030;
        21, 0.0150, 0.28, 0.0150, 0.032;
        23, 0.0630, 0.12, 0.0140, 0.023;
        25, 0.0300, 0.21, 0.0140, 0.031;
        28, 0.0220, 0.25, 0.0150, 0.033;
        30, 0.0660, 0.12, 0.0140, 0.022;
        40, 0.0260, 0.22, 0.0150, 0.029;
        42, 0.0089, 0.41, 0.0160, 0.043;
        64, 0.0120, 0.36, 0.0140, 0.047;
        65, 0.0048, 0.46, 0.0130, 0.073;
        79, 0.0350, 0.17, 0.0140, 0.026;
        94, 0.0280, 0.21, 0.0130, 0.034;
        100, 0.0130, 0.29, 0.0160, 0.033;
        101, 0.0079, 0.32, 0.0120, 0.074;
        139, 0.0240, 0.22, 0.0150, 0.028;
        147, 0.0470, 0.15, 0.0140, 0.026;
        173, 0.0390, 0.17, 0.0140, 0.029;
        175, 0.0057, 0.43, 0.0120, 0.073;
        200, 0.0110, 0.34, 0.0140, 0.053 
    ];

% Define solver options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    
% Plot the Pfizer vaccine data
    subplot(1, 2, 1); 
    hold on;
    
    for i = 1:size(pfizer_data, 1)
% Extract each patient's parameter from the Pfizer data
        r1 = pfizer_data(i, 2); 
        r2 = pfizer_data(i, 3);
        r3 = pfizer_data(i, 4); 
        r4 = pfizer_data(i, 5);
        
% Solve using ode45 
        [t, y] = ode45(@(t, y) system_odes(t, y, r1, r2, r3, r4, k1, k2, alpha, pfizer_days), tspan, y0, options);
        
 % Then plot V(t) 
        plot(t, y(:, 2), 'r', 'LineWidth', 1.0);
    end
    
    xlabel('Days Post-Vaccination'); 
    ylabel('Transfected Cells Proxy, V(t) (mL)');
    title('Pfizer Cohort Expression (Doses: Day 0, 21)');
    grid on;
    hold off;

% Plot the moderna data
    subplot(1, 2, 2); 
    hold on;
    
    for i = 1:size(moderna_data, 1)
% Extract each patient's parameter
        r1 = moderna_data(i, 2); 
        r2 = moderna_data(i, 3);
        r3 = moderna_data(i, 4); 
        r4 = moderna_data(i, 5);
        
% Solve using ode45
        [t, y] = ode45(@(t, y) system_odes(t, y, r1, r2, r3, r4, k1, k2, alpha, moderna_days), tspan, y0, options);
        
 % Then plot V(t)
        plot(t, y(:, 2), 'b', 'LineWidth', 1.0); 
    end
    
    xlabel('Days Post-Vaccination'); 
    ylabel('Transfected Cells Proxy, V(t) (mL)');
    title('Moderna Cohort Expression (Doses: Day 0, 28)');
    grid on;
    hold off;
    
% Main title for the entire figure
    sgtitle('Comparative Transfected Cell Kinetic Profiles V(t)', 'FontSize', 12, 'FontWeight', 'bold');
end

% non linear ODEs 
function dydt = system_odes(t, y, r1, r2, r3, r4, k1, k2, alpha, vaccine_days)
% State variables
    A = y(1); 
    V = y(2);
    
    % Determine if a vaccine pulse is currently active
    u_t = 0;
    for j = 1:length(vaccine_days)
        day = vaccine_days(j);
        if t >= day && t <= (day + 1.0) 
            u_t = 1.0; 
            break; 
        end
    end
    
    % The differential equations
    dAdt = (r1 * V) + (r2 * A * V) + A * (r3 - (r4 * A));
    dVdt = (alpha * u_t) - ((k1 * V) / (k2 + V));
    
    dydt = [dAdt; dVdt];
end