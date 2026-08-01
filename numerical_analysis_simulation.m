function numerical_analysis_simulation()

    % MASTER APP INTERFACE INITIALIZATION
    fig = uifigure('Name', 'SARS-CoV-2 Vaccine Kinetics Simulator (RK4)', ...
                   'Position', [50, 50, 1250, 750], 'Color', '#F4F5F7');

    % Create structural sub-panels for clean organization
    controlPanelL = uipanel(fig, 'Title', 'System Kinematics Controls', ...
                            'Position', [20, 20, 320, 710], 'FontWeight', 'bold', 'BackgroundColor', 'w');
    controlPanelR = uipanel(fig, 'Title', 'Patient Specific Vectors (r)', ...
                            'Position', [360, 20, 320, 710], 'FontWeight', 'bold', 'BackgroundColor', 'w');
    plotPanel     = uipanel(fig, 'Title', 'Dynamic Metric Visualization Engine', ...
                            'Position', [700, 20, 530, 710], 'FontWeight', 'bold', 'BackgroundColor', 'w');

    % Allocating Axis Objects Inside the Plot Panel 
    ax1 = uiaxes(plotPanel, 'Position', [40, 480, 460, 200]);
    ax2 = uiaxes(plotPanel, 'Position', [40, 260, 460, 160]);
    ax3 = uiaxes(plotPanel, 'Position', [40, 40, 460, 180]);

    % INITIAL ATTRIBUTE DECLARATION & INTERACTION STORAGE
    
    % Predefined baseline parameter states
    state.k1 = 10.0;   state.k2 = 50.0;   state.alpha = 1.0;
    state.a0 = 0.1;    state.v0 = 0.0;
    state.r1 = 0.0062; state.r2 = 0.40;   state.r3 = 0.0072; state.r4 = 0.027;

    % Core timeline setup 
    state.t = 0:0.1:45;
    
    % Building interactive sliders sequentially inside left panel
    sld_k1    = create_slider(controlPanelL, [30, 600, 260, 45], 'k1 Clearance Rate', [1, 50], state.k1);
    sld_k2    = create_slider(controlPanelL, [30, 500, 260, 45], 'k2 Half-Sat Constant', [5, 150], state.k2);
    sld_alpha = create_slider(controlPanelL, [30, 400, 260, 45], 'alpha Transcription Multiplier', [0.1, 5.0], state.alpha);
    sld_a0    = create_slider(controlPanelL, [30, 250, 260, 45], 'a0 Starting Antibodies', [0.0, 2.0], state.a0);
    sld_v0    = create_slider(controlPanelL, [30, 140, 260, 45], 'v0 Starting Transfected Cells', [0.0, 5.0], state.v0);

    % Building interactive sliders sequentially inside right panel
    sld_r1    = create_slider(controlPanelR, [30, 600, 260, 45], 'r1 Antibody Production Rate', [0.001, 0.05], state.r1);
    sld_r2    = create_slider(controlPanelR, [30, 450, 260, 45], 'r2 Neut. Elimination Factor', [0.01, 1.5], state.r2);
    sld_r3    = create_slider(controlPanelR, [30, 300, 260, 45], 'r3 Homeostatic Growth Rate', [0.001, 0.05], state.r3);
    sld_r4    = create_slider(controlPanelR, [30, 150, 260, 45], 'r4 Degradation Saturation', [0.005, 0.1], state.r4);

    
    % MAIN NUMERICAL INTEGRATION RUNNER (COMPUTATIONAL ENGINE)
    
    function run_simulation()
        t_vec = state.t;
        num_steps = length(t_vec);
        dt = t_vec(2) - t_vec(1);
        
        % Memory Preallocation
        a_rk4 = zeros(1, num_steps);
        v_rk4 = zeros(1, num_steps);
        u_vec = zeros(1, num_steps);
        
        a_rk4(1) = state.a0;
        v_rk4(1) = state.v0;
        
        % Inline Vectorized Step Pulse Condition
        evaluate_u = @(tau) double((tau >= 0.0 && tau <= 1.0) || (tau >= 21.0 && tau <= 22.0));
        u_vec(1) = evaluate_u(t_vec(1));
        
        % 4th Order Runge-Kutta Loop Mechanics
        for i = 1:(num_steps - 1)
            t_now = t_vec(i); a_now = a_rk4(i); v_now = v_rk4(i);
            
            % Stage 1
            u1 = evaluate_u(t_now);
            dadt1 = state.r1 * v_now + state.r2 * a_now * v_now + a_now * (state.r3 - state.r4 * a_now);
            dvdt1 = state.alpha * u1 - (state.k1 * v_now) / (state.k2 + v_now);
            
            % Stage 2
            t_mid = t_now + dt/2;
            a_mid1 = a_now + dadt1 * (dt / 2); v_mid1 = v_now + dvdt1 * (dt / 2);
            u2 = evaluate_u(t_mid);
            dadt2 = state.r1 * v_mid1 + state.r2 * a_mid1 * v_mid1 + a_mid1 * (state.r3 - state.r4 * a_mid1);
            dvdt2 = state.alpha * u2 - (k1_calc(v_mid1) / (state.k2 + v_mid1));
            
            % Stage 3
            a_mid2 = a_now + dadt2 * (dt / 2); v_mid2 = v_now + dvdt2 * (dt / 2);
            u3 = evaluate_u(t_mid);
            dadt3 = state.r1 * v_mid2 + state.r2 * a_mid2 * v_mid2 + a_mid2 * (state.r3 - state.r4 * a_mid2);
            dvdt3 = state.alpha * u3 - (k1_calc(v_mid2) / (state.k2 + v_mid2));
            
            % Stage 4
            t_end_step = t_now + dt;
            a_end = a_now + dadt3 * dt; v_end = v_now + dvdt3 * dt;
            u4 = evaluate_u(t_end_step);
            dadt4 = state.r1 * v_end + state.r2 * a_end * v_end + a_end * (state.r3 - state.r4 * a_end);
            dvdt4 = state.alpha * u4 - (k1_calc(v_end) / (state.k2 + v_end));
            
            % Final Combination Step
            a_rk4(i + 1) = a_now + (dt / 6) * (dadt1 + 2*dadt2 + 2*dadt3 + dadt4);
            v_rk4(i + 1) = v_now + (dt / 6) * (dvdt1 + 2*dvdt2 + 2*dvdt3 + dvdt4);
            u_vec(i + 1) = evaluate_u(t_vec(i + 1));
        end
        
        state.a_rk4 = a_rk4;
        state.v_rk4 = v_rk4;
        state.u_vec = u_vec;
    end

    function out = k1_calc(v_val)
        out = state.k1 * v_val;
    end

    
    % STATIC LAYOUT RENDERING & HANDLE CREATION
    
    run_simulation();
    
    % Initializing structural components for Axes 1 (Dual Y-Axis Emulation)
    yyaxis(ax1, 'left');
    h_a = plot(ax1, state.t, state.a_rk4, 'b-', 'LineWidth', 2);
    ylabel(ax1, 'Antibodies a(t)', 'Color', 'b', 'FontWeight', 'bold');
    ax1.YColor = 'b';
    grid(ax1, 'on');
    
    % Realtime Peak Tracking Handles
    [peak_val, peak_idx] = max(state.a_rk4);
    h_pk = line(ax1, state.t(peak_idx), peak_val, 'Marker', 'o', ...
                'MarkerFaceColor', 'r', 'Color', 'r', 'MarkerSize', 6);
    
    yyaxis(ax1, 'right');
    h_v = plot(ax1, state.t, state.v_rk4, 'r--', 'LineWidth', 1.8);
    ylabel(ax1, 'Transfected Cells v(t)', 'Color', 'r', 'FontWeight', 'bold');
    ax1.YColor = 'r';
    title(ax1, 'Immune Engine Tracking Profile');
    
    % Initializing Axes 2 (Dose Graph Area)
    h_u = area(ax2, state.t, state.u_vec, 'FaceColor', [0.2 0.7 0.3], ...
               'EdgeColor', [0.1 0.5 0.2], 'FaceAlpha', 0.4);
    grid(ax2, 'on');
    ylabel(ax2, 'Dose Profile u(t)');
    xlabel(ax2, 'Time (Days)');
    ylim(ax2, [-0.1 1.2]);
    
    % Initializing Axes 3 (Phase Space Profile)
    h_ps   = plot(ax3, state.v_rk4, state.a_rk4, 'm-', 'LineWidth', 2);
    grid(ax3, 'on'); hold(ax3, 'on');
    h_ps_s = plot(ax3, state.v_rk4(1), state.a_rk4(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    h_ps_e = plot(ax3, state.v_rk4(end), state.a_rk4(end), 'rx', 'LineWidth', 2, 'MarkerSize', 9);
    xlabel(ax3, 'Transfected Cells v(t)');
    ylabel(ax3, 'Antibody Level a(t)');
    title(ax3, 'Phase Space State Boundary Map');

    
    % REALTIME VALUE CHANGE DISPATCHER CALLBACKS
    function respond_callback(~, event, parameter_id)
        % Bind changing UI scalar values directly to engine structure
        state.(parameter_id) = event.Value;
        
        % Force update calculation run
        run_simulation();
        
        % Dispatch graphic structural transformations inside internal buffers
        h_a.YData = state.a_rk4;
        h_v.YData = state.v_rk4;
        h_u.YData = state.u_vec;
        
        [p_v, p_idx] = max(state.a_rk4);
        h_pk.XData = state.t(p_idx);
        h_pk.YData = p_v;
        
        h_ps.XData = state.v_rk4;  h_ps.YData = state.a_rk4;
        h_ps_s.XData = state.v_rk4(1); h_ps_s.YData = state.a_rk4(1);
        h_ps_e.XData = state.v_rk4(end); h_ps_e.YData = state.a_rk4(end);
    end

    % Helper design layout function to maintain clean encapsulation
    function sld = create_slider(parent, pos, display_label, limits, start_val)
        uilabel(parent, 'Position', [pos(1), pos(2)+25, pos(3), 20], ...
                'Text', display_label, 'FontWeight', 'bold');
        sld = uislider(parent, 'Position', pos, 'Limits', limits, 'Value', start_val);
        
        % Extract field character arrays using structural regex tokens
        clean_id = regexp(display_label, '^\w+', 'match', 'once');
        sld.ValueChangingFcn = @(src, event) respond_callback(src, event, clean_id);
    end
end
