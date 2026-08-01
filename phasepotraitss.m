%% 
% Phase Portrait - Antibodies vs. Transfected Cells
clear; clc; close all;
% Parameters & Initial Conditions
t1_span = [0 21];
t2_span = [21 180]; % Extended to 180 days to show decay down y-axis to steady state
k1 = 10; k2 = 50; 
r1 = 0.019; r2 = 0.26; r3 = 0.014; r4 = 0.038;
% ODE System (Natural dynamics without continuous pulse)
sys_ode = @(t,y) [ ...
    r1*y(2) + r2*y(1)*y(2) + y(1)*(r3 - r4*y(1)); ...
    -(k1*y(2)) / (k2 + y(2) + eps) ...
];
options = odeset('MaxStep', 0.05); 
% Dose 1 (Day 0) - Instantaneous jump of V from 0 to 1 at A = 0.1
y0_dose1 = [0.1; 1]; 
[t1, y1] = ode45(sys_ode, t1_span, y0_dose1, options);
% Dose 2 (Day 21) - Instantaneous jump adding +1 to V at current antibody level
A_dose2_start = y1(end, 1);
V_dose2_start = y1(end, 2) + 1;
y0_dose2 = [A_dose2_start; V_dose2_start];
[t2, y2] = ode45(sys_ode, t2_span, y0_dose2, options);
% Combine data and insert horizontal jump coordinates
V_full = [0; y1(:,2); V_dose2_start; y2(:,2)];
A_full = [0.1; y1(:,1); A_dose2_start; y2(:,1)];
% Phase Portrait (V on X-axis, A on Y-axis)
figure(10);
set(gcf, 'Color', 'w');
plot(V_full, A_full, 'k-', 'LineWidth', 1.5);
title('phase portrait: antibodies vs. transfected cells', 'FontWeight', 'bold');
xlabel('transfected cells proxy v(t)');
ylabel('antibody level a(t)');
grid on; 
grid minor;
xlim([0, 1.2]);
ylim([0, 3.5]);