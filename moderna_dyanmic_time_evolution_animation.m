
close all force; 
clear; clc;

%  Parameters 
k1 = 10.0; 
k2 = 50.0; 
alpha = 1.0;
a0 = 0.1;  
v0 = 0.0;
r = [0.02, 0.2, 0.01, 0.03]; % Baseline parameters [r1, r2, r3, r4]
schedule = [0.0, 28.0];      % MODERNA SCHEDULE: Dose 1 at Day 0, Dose 2 at Day 28

% Simulation Time Array 
dt = 0.05; 
t = 0:dt:45;
num_steps = length(t);

% Run High-Precision RK4 Forward Integration 
a = zeros(1, num_steps); v = zeros(1, num_steps);
a(1) = a0; v(1) = v0;

for i = 1:(num_steps - 1)
    u_t = 0.0;
    if (t(i) >= schedule(1) && t(i) <= schedule(1)+1.0) || (t(i) >= schedule(2) && t(i) <= schedule(2)+1.0)
        u_t = 1.0;
    end
    
    % Stage 1
    da1 = r(1)*v(i) + r(2)*a(i)*v(i) + a(i)*(r(3) - r(4)*a(i));
    dv1 = alpha*u_t - (k1*v(i))/(k2 + v(i));
    
    % Stage 2
    am1 = a(i) + da1*(dt/2); vm1 = v(i) + dv1*(dt/2);
    da2 = r(1)*vm1 + r(2)*am1*vm1 + am1*(r(3) - r(4)*am1);
    dv2 = alpha*u_t - (k1*vm1)/(k2 + vm1);
    
    % Stage 3
    am2 = a(i) + da2*(dt/2); vm2 = v(i) + dv2*(dt/2);
    da3 = r(1)*vm2 + r(2)*am2*vm2 + am2*(r(3) - r(4)*am2);
    dv3 = alpha*u_t - (k1*vm2)/(k2 + vm2);
    
    % Stage 4
    ae = a(i) + da3*dt; ve = v(i) + dv3*dt;
    da4 = r(1)*ve + r(2)*ae*ve + ae*(r(3) - r(4)*ae);
    dv4 = alpha*u_t - (k1*ve)/(k2 + ve);
    
    % Final Weighted Step
    a(i+1) = a(i) + (dt/6)*(da1 + 2*da2 + 2*da3 + da4);
    v(i+1) = v(i) + (dt/6)*(dv1 + 2*dv2 + 2*dv3 + dv4);
end


video_filename = 'Moderna_System_Dynamics.mp4';
v_out = VideoWriter(video_filename, 'MPEG-4');
v_out.FrameRate = 15;  
v_out.Quality = 95;    

try
    open(v_out);

    % establish figure canvas properties
    fig = figure(1);
    set(fig, 'Position', [100, 100, 960, 540], 'Color', 'w'); 

    ax = axes('Parent', fig);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'Layer', 'top'); 
    xlim(ax, [0, 45]); ylim(ax, [-0.05, 1.25]);
    xlabel(ax, 'Timeline (Days)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Relative Concentration Level', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Updated Explicit Title
    title(ax, 'Moderna Cohort Trajectory Simulation: a(t) vs v(t) (Booster Day 28)', 'FontSize', 14, 'FontWeight', 'bold');

    % Duller, highly-transparent background panels for visibility
    mute_color = [0.94, 0.94, 0.95]; 
    patch([0 1 1 0], [-0.05 -0.05 1.25 1.25], mute_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Dose Active Windows');
    patch([28 29 29 28], [-0.05 -0.05 1.25 1.25], mute_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

    % Data Line Handles
    h_antigen  = plot(ax, NaN, NaN, 'Color', [0.95, 0.45, 0.1], 'LineWidth', 2.5, 'DisplayName', 'Antigen v(t)');
    h_antibody = plot(ax, NaN, NaN, 'Color', [0.1, 0.55, 0.85], 'LineWidth', 2.5, 'DisplayName', 'Antibody a(t)');
    h_cursor   = plot(ax, [NaN, NaN], [-0.05, 1.25], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    legend(ax, 'Location', 'northeast', 'FontSize', 11);

    % Animation compilation loop
    frame_skip = 2; 

    for idx = 1:frame_skip:num_steps
        if ~ishandle(fig), break; end
        
        set(h_antigen, 'XData', t(1:idx), 'YData', v(1:idx));
        set(h_antibody, 'XData', t(1:idx), 'YData', a(1:idx));
        set(h_cursor, 'XData', [t(idx), t(idx)]);
        
        drawnow;
        
        frame = getframe(fig);
        writeVideo(v_out, frame);
    end

    close(v_out);
    fprintf('\n>>> Moderna video compilation completed successfully! Saved as: %s\n', video_filename);

catch ME
    close(v_out);
    rethrow(ME);
end