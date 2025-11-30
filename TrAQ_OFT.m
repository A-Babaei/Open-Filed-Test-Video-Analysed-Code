% A.Babaeic
%% Robust TrAQ Open Field Analysis: Locomotion + Zones (Peripheral/Central)
% Requires: Ensure "Outut.mat" is in your current MATLAB directory.

arena_cm = 70; % real arena dimension (cm)
peripheral_width_cm = 10; % cm-wide peripheral border

%% AUTO-LOAD if needed
if ~exist('track', 'var')
    if ~exist('Out_Stim2', 'file')
        error('''Out_Stim2'' not found in current directory.');
    end
    load('Out_Stim2');
    fprintf('Loaded Out_Stim2\n');
end

%% DIAGNOSTIC CHECKS
if ~exist('track', 'var')
    error('''track'' not loaded. Run load(''Out_Stimulation.mat'');');
end
if ~isstruct(track)
    error('''track'' is not a struct (class: %s).', class(track));
end
fprintf('track is a struct with fields: %s\n', strjoin(fieldnames(track), ', '));

%% Set default font for figures
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontWeight', 'bold');

%% 1. Determine frame rate and time base - IMPROVED
if isfield(track, 'Time') && ~isempty(track.Time)
    Time = track.Time(:)';
    if length(Time) > 1
        dt = mean(diff(Time));
        fps = 1 / dt;
    else
        dt = 1/30;
        fps = 30;
    end
    nFrames = numel(Time);
    time_source = 'track.Time';
elseif isfield(track, 'FrameRate')
    fps = track.FrameRate;
    dt = 1 / fps;
    nFrames = size(track.Centroid, 1);
    if nFrames < 100
        nFrames = size(track.Centroid, 2);
    end
    time_source = 'track.FrameRate';
else
    fps = 30;
    dt = 1/fps;
    nFrames = max(size(track.Centroid));
    time_source = 'default(30)';
end
total_time_s = nFrames * dt;
fprintf('Detected nFrames: %d, fps: %.2f, total time: %.2f s\n', nFrames, fps, total_time_s);

%% 2. Load positions for movement and zones
C = []; position_source = '';
if isfield(track, 'Centroid')
    C = double(track.Centroid);
    position_source = 'track.Centroid';
    fprintf('Using %s for positions (size: %dx%d, nFrames=%d)\n', position_source, size(C,1), size(C,2), nFrames);
elseif isfield(track, 'Head')
    C = double(track.Head);
    position_source = 'track.Head';
    fprintf('Using %s for positions (size: %dx%d, nFrames=%d)\n', position_source, size(C,1), size(C,2), nFrames);
else
    C_fields = {'Tail', 'Body_Centroid', 'Position', 'HeadPosition', 'NosePosition', 'BodyPosition'};
    for i = 1:length(C_fields)
        if isfield(track, C_fields{i})
            temp_C = double(track.(C_fields{i}));
            if ((size(temp_C, 1) == nFrames && size(temp_C, 2) == 2) || (size(temp_C, 2) == nFrames && size(temp_C, 1) == 2))
                C = temp_C;
                position_source = ['track.', C_fields{i}];
                fprintf('Using %s for positions\n', position_source);
                break;
            end
        end
    end
end

if isempty(C)
    error('No suitable position data found.');
end

%% Data validation and cleaning
fprintf('\n--- Data Validation ---\n');
fprintf('Position data range: X [%.1f to %.1f], Y [%.1f to %.1f] px\n', ...
    min(C(:,1)), max(C(:,1)), min(C(:,2)), max(C(:,2)));

% Remove or interpolate NaN values
nan_mask = any(isnan(C), 2);
if sum(nan_mask) > 0
    fprintf('Warning: %d frames have NaN positions - interpolating\n', sum(nan_mask));
    for col = 1:2
        C(:,col) = fillmissing(C(:,col), 'linear');
    end
end

%% Ensure C is N×2 for zone analysis
if size(C, 1) == 2
    C = C';
end
C = C(1:nFrames, :);

%% Compute distance per frame from positions - IMPROVED
if size(C,1) == 2
    % 2×N matrix (x;y)
    dx = diff(C(1,:));
    dy = diff(C(2,:));
    d_px = [0, sqrt(dx.^2 + dy.^2)];
else
    % Nx2 matrix [x,y]
    dx = diff(C(:,1));
    dy = diff(C(:,2));
    d_px = [0; sqrt(dx.^2 + dy.^2)]';
end

% Filter out unrealistic jumps (tracking errors)
max_reasonable_jump_px = 30; % Adjust based on your data
large_jumps = d_px > max_reasonable_jump_px;
if sum(large_jumps) > 0
    fprintf('Filtering %d unrealistic position jumps (>%.1f px)\n', sum(large_jumps), max_reasonable_jump_px);
    d_px(large_jumps) = 0;
end

d_px = d_px(1:nFrames);

%% 3. Pixel → cm conversion
arena_available = false;
if exist('arena', 'var') && ~isempty(arena)
    arena_available = true;
    min_x = min(arena(:,1)); max_x = max(arena(:,1));
    min_y = min(arena(:,2)); max_y = max(arena(:,2));
    width_px = max_x - min_x;
    height_px = max_y - min_y;
    arena_diameter_px = mean([width_px, height_px]);
    px_to_cm = arena_cm / arena_diameter_px;
    fprintf('Arena: %.1f cm = %.1f px → 1 px = %.5f cm\n', arena_cm, arena_diameter_px, px_to_cm);
else
    warning('Arena not loaded; using estimated conversion');
    px_to_cm = 70/500;
    fprintf('Using estimated conversion: 1 px = %.5f cm\n', px_to_cm);
    min_x = min(C(:,1)); max_x = max(C(:,1));
    min_y = min(C(:,2)); max_y = max(C(:,2));
    width_px = max_x - min_x;
    height_px = max_y - min_y;
end

%% 4. Convert distance and compute speed - WITH REALISTIC LIMITS
d_cm = d_px * px_to_cm;
speed_cm_s = d_cm / dt;

% Debug speed calculation
fprintf('\n--- Speed Validation ---\n');
fprintf('Max distance between frames: %.3f cm\n', max(d_cm));
fprintf('Frame duration (dt): %.6f seconds\n', dt);
fprintf('Theoretical max speed if 10px/frame: %.3f cm/s\n', (10 * px_to_cm) / dt);
fprintf('Actual max speed before filtering: %.3f cm/s\n', max(speed_cm_s));

% Apply realistic speed limits for rodent locomotion
max_realistic_speed = 60; % cm/s - maximum realistic speed for rodents
speed_cm_s = min(speed_cm_s, max_realistic_speed);

% Apply mild smoothing to reduce noise
speed_cm_s = smoothdata(speed_cm_s, 'movmean', 3);

fprintf('Final speed stats - max: %.2f, 95th percentile: %.2f, mean: %.2f cm/s\n', ...
    max(speed_cm_s), prctile(speed_cm_s, 95), mean(speed_cm_s));

total_relocation_cm = sum(d_cm);

%% 5. Movement threshold (define locomotion)
movement_threshold = 2.0; % cm/s
moving_idx = speed_cm_s > movement_threshold;

%% 6. Compute time in locomotion and immobility
time_moving_s = sum(moving_idx) * dt;
time_immobile_s = total_time_s - time_moving_s;
percent_moving = 100 * sum(moving_idx) / nFrames;

%% 7. Speed summaries
mean_speed_all = mean(speed_cm_s);
mean_speed_moving = mean(speed_cm_s(moving_idx));

%% 8. ZONE ANALYSIS: Peripheral vs. Central
peripheral_width_px = peripheral_width_cm / px_to_cm;
central_mask = (C(:,1) > min_x + peripheral_width_px) & (C(:,1) < max_x - peripheral_width_px) & ...
               (C(:,2) > min_y + peripheral_width_px) & (C(:,2) < max_y - peripheral_width_px);
peripheral_mask = ~central_mask;

time_central_s = sum(central_mask) * dt;
time_peripheral_s = sum(peripheral_mask) * dt;
percent_central = 100 * time_central_s / total_time_s;
percent_peripheral = 100 * time_peripheral_s / total_time_s;

mean_speed_central = mean(speed_cm_s(central_mask));
mean_speed_peripheral = mean(speed_cm_s(peripheral_mask));
total_dist_central_cm = sum(d_cm(central_mask));
total_dist_peripheral_cm = sum(d_cm(peripheral_mask));

%% 9. Display results
fprintf('\n--- Open Field Locomotion Summary ---\n');
fprintf('Frame rate: %.2f fps (source: %s)\n', fps, time_source);
fprintf('Total frames: %d | Total time: %.2f s\n', nFrames, total_time_s);
fprintf('Total distance moved: %.2f cm\n', total_relocation_cm);
fprintf('Total time moving: %.2f s (%.1f%%)\n', time_moving_s, percent_moving);
fprintf('Mean speed (all frames): %.3f cm/s\n', mean_speed_all);
fprintf('Mean speed (moving only): %.3f cm/s\n', mean_speed_moving);
fprintf('\n--- Zone Summary ---\n');
fprintf('Time in Central: %.2f s (%.1f%%)\n', time_central_s, percent_central);
fprintf('Time in Peripheral: %.2f s (%.1f%%)\n', time_peripheral_s, percent_peripheral);
fprintf('Mean speed in Central: %.3f cm/s\n', mean_speed_central);
fprintf('Mean speed in Peripheral: %.3f cm/s\n', mean_speed_peripheral);

%% 10. Visualization
t = (0:nFrames-1) * dt;

% Figure 1: Speed plots
fig1 = figure('Name','Speed Plot','NumberTitle','off');
subplot(2,1,1);
plot(t, speed_cm_s, 'b-', 'LineWidth', 1); hold on;
yline(movement_threshold, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Threshold');
xlabel('Time [s]'); ylabel('Speed [cm/s]');
title('Speed Plot'); grid on;
ylim([0 60]); % Set realistic y-axis limits
legend('Speed', 'Threshold', 'Location', 'best');

subplot(2,1,2);
histogram(speed_cm_s, 50, 'Normalization', 'probability');
xlabel('Speed [cm/s]'); ylabel('Probability [%]');
title('Speed Probability'); grid on;
xline(movement_threshold, 'r--', 'Threshold');
xlim([0 60]); % Set realistic x-axis limits

% Figure 2: Zone analysis
fig2 = figure('Name','Zone Occupancy','NumberTitle','off');

% Fixed histogram labels
subplot(2,2,1); 
histogram(speed_cm_s(central_mask), 50); 
title('Speed in Central');
xlabel('Speed (cm/s)'); ylabel('Frequency');
xline(movement_threshold, 'r--');
xlim([0 60]);

subplot(2,2,2); 
histogram(speed_cm_s(peripheral_mask), 50); 
title('Speed in Peripheral');
xlabel('Speed (cm/s)'); ylabel('Frequency');
xline(movement_threshold, 'r--');
xlim([0 60]);

subplot(2,2,3);
stairs(t, cumsum(central_mask)*dt, 'g-', 'LineWidth', 2); hold on;
stairs(t, cumsum(peripheral_mask)*dt, 'm-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Cumulative Time (s)');
title('Zone Time Accumulation');
legend('Central', 'Peripheral', 'Location', 'best');

% Trajectory plot in PIXELS
subplot(2,2,4);
x_px = C(:,1);
y_px = C(:,2);

plot(x_px, y_px, 'k-', 'LineWidth', 1, 'DisplayName', 'Animal Path');
hold on;

% Peripheral zone
peripheral_x_px = [min_x, max_x, max_x, min_x, min_x];
peripheral_y_px = [min_y, min_y, max_y, max_y, min_y];
patch(peripheral_x_px, peripheral_y_px, 'r', 'FaceAlpha', 0, ...
      'EdgeColor', 'r', 'LineWidth', 3, 'DisplayName', 'Peripheral Zone');

% Central zone
central_start_x_px = min_x + peripheral_width_px;
central_start_y_px = min_y + peripheral_width_px;
central_end_x_px = max_x - peripheral_width_px;
central_end_y_px = max_y - peripheral_width_px;

if (central_end_x_px > central_start_x_px) && (central_end_y_px > central_start_y_px)
    central_x_px = [central_start_x_px, central_end_x_px, central_end_x_px, central_start_x_px, central_start_x_px];
    central_y_px = [central_start_y_px, central_start_y_px, central_end_y_px, central_end_y_px, central_start_y_px];
    patch(central_x_px, central_y_px, 'g', 'FaceAlpha', 0, ...
          'EdgeColor', 'g', 'LineWidth', 2, 'DisplayName', 'Central Zone');
end

xlim([min_x, max_x]);
ylim([min_y, max_y]);
xlabel('X [pix]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y [pix]', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Trajectory (Peripheral = %.1f cm)', peripheral_width_cm));
grid on;
axis equal;
legend('Location', 'best', 'FontSize', 10);

%% Save results
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
baseFilename = sprintf('TrAQ_Results_%s', timestamp);
saveas(fig1, sprintf('%s_speed_plot.fig', baseFilename));
saveas(fig1, sprintf('%s_speed_plot.bmp', baseFilename));
saveas(fig2, sprintf('%s_zones_diagnostics.fig', baseFilename));
saveas(fig2, sprintf('%s_zones_diagnostics.bmp', baseFilename));

fprintf('\nResults saved with timestamp: %s\n', timestamp);
fprintf('Speed values now capped at realistic maximum (60 cm/s)\n');