% =========================================================================
% MATLAB ANIMATION EXPORTER: COUPLED IMMUNE DYNAMICS (MP4 GENERATOR)
% =========================================================================
% Description: Integrates the non-linear a(t)/v(t) system and compiles a
% frame-by-frame animation exported directly as a high-fidelity MP4 video.
% =========================================================================

clear; clc; close all;

% --- Core Physical Parameters ---
k1 = 10.0; k2 = 50.0; alpha = 1.0;
a0 = 0.1;  v0 = 0.0;
r1 = 0.02; r2 = 0.2; r3 = 0.01; r4 = 0.03; % Individual parameters
schedule = [0.0, 21.0];      % Boost schedules (Days)

% --- Simulation Time Array ---
dt = 0.05; % Fine step size for smooth animation frames
t = 0:dt:45;
num_steps = length(t);

% =========================================================================
% CORRECTED RUNGE-KUTTA 4TH ORDER (RK4) LOOP IMPLEMENTATION
% =========================================================================

% Anonymous helper function to evaluate the pulse function dynamically at any given time 'tau'
evaluate_u = @(tau) double((tau >= 0.0 && tau <= 1.0) || (tau >= 21.0 && tau <= 22.0));

% Reset storage and initial states
a = zeros(1, num_steps);
v = zeros(1, num_steps);
a(1) = a0;
v(1) = v0;

for i = 1:(num_steps - 1)
    t_now = t(i);
    a_now = a(i);
    v_now = v(i);
    
    % --- STAGE 1: Evaluated at the start of the interval (t) ---
    u1 = evaluate_u(t_now);
    dadt1 = r1 * v_now + r2 * a_now * v_now + a_now * (r3 - r4 * a_now);
    dvdt1 = alpha * u1 - (k1 * v_now) / (k2 + v_now);
    
    % --- STAGE 2: Evaluated at the first midpoint (t + dt/2) ---
    t_mid = t_now + dt/2;
    a_mid1 = a_now + dadt1 * (dt / 2);
    v_mid1 = v_now + dvdt1 * (dt / 2);
    
    u2 = evaluate_u(t_mid); % Recalculate pulse at midpoint!
    dadt2 = r1 * v_mid1 + r2 * a_mid1 * v_mid1 + a_mid1 * (r3 - r4 * a_mid1);
    dvdt2 = alpha * u2 - (k1 * v_mid1) / (k2 + v_mid1);
    
    % --- STAGE 3: Evaluated at the second midpoint (t + dt/2) ---
    a_mid2 = a_now + dadt2 * (dt / 2);
    v_mid2 = v_now + dvdt2 * (dt / 2);
    
    u3 = evaluate_u(t_mid); % Recalculate pulse at midpoint!
    dadt3 = r1 * v_mid2 + r2 * a_mid2 * v_mid2 + a_mid2 * (r3 - r4 * a_mid2);
    dvdt3 = alpha * u3 - (k1 * v_mid2) / (k2 + v_mid2);
    
    % --- STAGE 4: Evaluated at the end of the interval (t + dt) ---
    t_end = t_now + dt;
    a_end = a_now + dadt3 * dt;
    v_end = v_now + dvdt3 * dt;
    
    u4 = evaluate_u(t_end); % Recalculate pulse at endpoint!
    dadt4 = r1 * v_end + r2 * a_end * v_end + a_end * (r3 - r4 * a_end);
    dvdt4 = alpha * u4 - (k1 * v_end) / (k2 + v_end);
    
    % --- FINAL WEIGHTED COMBINATION ---
    a(i + 1) = a_now + (dt / 6) * (dadt1 + 2*dadt2 + 2*dadt3 + dadt4);
    v(i + 1) = v_now + (dt / 6) * (dvdt1 + 2*dvdt2 + 2*dvdt3 + dvdt4);
end

% =========================================================================
% VIDEO WRITER INITIALIZATION
% =========================================================================
video_filename = 'Immune_System_Dynamics.mp4';
v_out = VideoWriter(video_filename, 'MPEG-4');
v_out.FrameRate = 15;  % REDUCED FROM 30 TO 15 (HALF SPEED)
v_out.Quality = 95;    % High video quality bit-rate
open(v_out);

% --- Establish Figure Canvas Properties ---
fig = figure(1);
set(fig, 'Position', [100, 100, 960, 540], 'Color', 'w'); % Standard 16:9 frame aspect ratio

% Setup axes
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [0, 45]);
ylim(ax, [-0.05, 1.25]);
xlabel(ax, 'Timeline (Days)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'Relative Concentration Level', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, 'Pfizer Dynamic Trajectory Simulation: a(t) vs v(t)', 'FontSize', 14, 'FontWeight', 'bold');

% Draw Static Dose Windows for structural reference
patch(ax, [0 1 1 0], [-0.05 -0.05 1.25 1.25], [0.9 0.9 0.9], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Dose Active Windows');
patch(ax, [21 22 22 21], [-0.05 -0.05 1.25 1.25], [0.9 0.9 0.9], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% Initialize dynamic line handles
h_antigen  = plot(ax, NaN, NaN, 'Color', [0.95, 0.45, 0.1], 'LineWidth', 2.5, 'DisplayName', 'Antigen v(t)');
h_antibody = plot(ax, NaN, NaN, 'Color', [0.1, 0.55, 0.85], 'LineWidth', 2.5, 'DisplayName', 'Antibody a(t)');
h_cursor   = plot(ax, [NaN, NaN], [-0.05, 1.25], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');

legend(ax, 'Location', 'northeast', 'FontSize', 11);

% =========================================================================
% ANIMATION COMPILATION LOOP
% =========================================================================
% Skip indices systematically to convert frame rates cleanly
frame_skip = 1;  % CHANGED FROM 2 TO 1 FOR SMOOTHER PLAYBACK WITH HALF FRAME RATE

% Force figure to render before capturing frames
drawnow;

for idx = 1:frame_skip:num_steps
    % Append data up to current frame index
    set(h_antigen, 'XData', t(1:idx), 'YData', v(1:idx));
    set(h_antibody, 'XData', t(1:idx), 'YData', a(1:idx));
    
    % Update the timeline marker line position
    set(h_cursor, 'XData', [t(idx), t(idx)]);
    
    drawnow;
    
    % Capture the frame from the current figure window
    frame = getframe(gcf);
    writeVideo(v_out, frame);
end

% Clean wrap-up closure
close(v_out);
fprintf('\n>>> Video compilation completed successfully! Saved as: %s\n', video_filename);