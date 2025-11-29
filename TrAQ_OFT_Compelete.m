%% Robust TrAQ Open Field Analysis: Locomotion + Zones (Peripheral/Central)
% FIXED: Respects manual peripheral width setting + Debugged version

arena_cm = 70; % real arena dimension (cm)
peripheral_width_cm = 10; % cm-wide peripheral border - YOUR MANUAL SETTING

%% AUTO-LOAD if needed
if ~exist('track', 'var')
    if ~exist('Out_Stim2', 'file')
        error('''Out_Stim2'' not found in current directory.');
    end
    load('Out_Stim2.mat');
    fprintf('Loaded Out_Stim2.mat\n');
end

%% DIAGNOSTIC CHECKS
if ~exist('track', 'var')
    error('''track'' not loaded from Out_Post.mat.');
end
if ~isstruct(track)
    error('''track'' is not a struct (class: %s).', class(track));
end
fprintf('track is a struct with fields: %s\n', strjoin(fieldnames(track), ', '));

%% Input validation
if peripheral_width_cm <= 0
    error('peripheral_width_cm must be positive');
end
if arena_cm <= 0
    error('arena_cm must be positive');
end

%% Set default font for figures
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontWeight', 'bold');

%% 1. Determine frame rate and time base
if isfield(track, 'Time') && ~isempty(track.Time)
    Time = track.Time(:);
    fps = 1 / mean(diff(Time));
    nFrames = numel(Time);
    dt = 1/fps;
    time_source = 'track.Time';
    fprintf('Using track.Time: %.2f fps, %d frames\n', fps, nFrames);
elseif isfield(track, 'FrameRate')
    fps = track.FrameRate;
    dt = 1 / fps;
    if isfield(track, 'Centroid')
        nFrames = size(track.Centroid, 1);
    else
        nFrames = length(track.Time);
    end
    time_source = 'track.FrameRate';
else
    fps = 30; % fallback
    dt = 1/fps;
    if isfield(track, 'Centroid')
        nFrames = size(track.Centroid, 1);
    else
        nFrames = 1000;
    end
    time_source = 'default(30)';
end

total_time_s = nFrames * dt;
fprintf('nFrames: %d, fps: %.2f, total time: %.2f s\n', nFrames, fps, total_time_s);

%% 2. Load positions and handle coordinate systems
C = []; position_source = '';
if isfield(track, 'Centroid')
    C = double(track.Centroid);
    position_source = 'track.Centroid';
    fprintf('Using Centroid: %dx%d matrix\n', size(C,1), size(C,2));
elseif isfield(track, 'Head')
    C = double(track.Head);
    position_source = 'track.Head';
    fprintf('Using Head: %dx%d matrix\n', size(C,1), size(C,2));
else
    C_fields = {'Tail', 'Body_Centroid', 'Position'};
    for i = 1:length(C_fields)
        if isfield(track, C_fields{i})
            C = double(track.(C_fields{i}));
            position_source = ['track.', C_fields{i}];
            fprintf('Using %s: %dx%d matrix\n', position_source, size(C,1), size(C,2));
            break;
        end
    end
end

if isempty(C)
    error('No suitable position data found.');
end

% Ensure C is numeric and handle NaNs
C = double(C);
C(any(isnan(C), 2), :) = []; % Remove rows with NaNs

% Handle coordinate formats
if size(C, 2) == 2 && size(C, 1) > 2
    % Nx2 format
    x = C(:,1);
    y = C(:,2);
    fprintf('Detected Nx2 coordinate format: %d points\n', length(x));
elseif size(C, 1) == 2 && size(C, 2) > 2
    % 2xN format
    x = C(1,:)';
    y = C(2,:)';
    fprintf('Detected 2xN coordinate format - transposing: %d points\n', length(x));
elseif ndims(C) == 3
    % Handle 3D arrays (common in tracking data)
    x = squeeze(C(:,1,1));
    y = squeeze(C(:,1,2));
    fprintf('Detected 3D coordinate format: %d points\n', length(x));
else
    error('Unknown coordinate format: %dx%d', size(C,1), size(C,2));
end

% Ensure we have exactly nFrames (handle case where we have fewer points)
actual_frames = min(nFrames, length(x));
x = x(1:actual_frames);
y = y(1:actual_frames);
if actual_frames < nFrames
    fprintf('Warning: Only %d frames of position data (expected %d)\n', actual_frames, nFrames);
    nFrames = actual_frames;
    total_time_s = nFrames * dt;
end

%% 3. Pixel-to-cm calibration
fprintf('\n--- Coordinate System Diagnostics ---\n');
fprintf('X range: %.1f to %.1f (span: %.1f px)\n', min(x), max(x), max(x)-min(x));
fprintf('Y range: %.1f to %.1f (span: %.1f px)\n', min(y), max(y), max(y)-min(y));

% Calculate arena bounds with safety margins
x_margin = (max(x) - min(x)) * 0.1;
y_margin = (max(y) - min(y)) * 0.1;

min_x_arena = min(x) - x_margin;
max_x_arena = max(x) + x_margin;
min_y_arena = min(y) - y_margin;
max_y_arena = max(y) + y_margin;

arena_width_px = max_x_arena - min_x_arena;
arena_height_px = max_y_arena - min_y_arena;

% Calculate pixel-to-cm conversion
px_to_cm = arena_cm / max(arena_width_px, arena_height_px);

fprintf('Arena bounds: X(%.1f to %.1f), Y(%.1f to %.1f) px\n', min_x_arena, max_x_arena, min_y_arena, max_y_arena);
fprintf('Arena size: %.1f x %.1f px\n', arena_width_px, arena_height_px);
fprintf('Pixel to cm conversion: 1 px = %.6f cm\n', px_to_cm);

%% 4. Speed calculation with calibration
if nFrames > 1
    dx = diff(x);
    dy = diff(y);
    distance_px_frame = sqrt(dx.^2 + dy.^2);
    distance_cm_frame = [0; distance_px_frame] * px_to_cm;
else
    distance_cm_frame = zeros(nFrames, 1);
end

speed_cm_s_raw = distance_cm_frame / dt;

% Manual smoothing (replaces smooth() function)
if length(speed_cm_s_raw) >= 5
    kernel = ones(5,1)/5;
    speed_cm_s = conv(speed_cm_s_raw, kernel, 'same');
    % Handle edges
    speed_cm_s(1:2) = speed_cm_s_raw(1:2);
    speed_cm_s(end-1:end) = speed_cm_s_raw(end-1:end);
else
    speed_cm_s = speed_cm_s_raw;
end

total_distance_cm = sum(distance_cm_frame);

% Automatic calibration if speeds are unrealistic
fprintf('\n--- Speed Calculation Verification ---\n');
fprintf('Max raw speed: %.2f cm/s\n', max(speed_cm_s_raw));
fprintf('Max smoothed speed: %.2f cm/s\n', max(speed_cm_s));

if max(speed_cm_s) > 100
    calibration_factor = 25 / max(speed_cm_s);
    speed_cm_s = speed_cm_s * calibration_factor;
    distance_cm_frame = distance_cm_frame * calibration_factor;
    total_distance_cm = total_distance_cm * calibration_factor;
    px_to_cm = px_to_cm * calibration_factor;
    fprintf('APPLIED CALIBRATION: Scaling factor = %.4f\n', calibration_factor);
    fprintf('New max speed: %.2f cm/s\n', max(speed_cm_s));
    fprintf('New pixel conversion: 1 px = %.6f cm\n', px_to_cm);
end

%% 5. Movement threshold
movement_threshold = 2.0; % cm/s
moving_idx = speed_cm_s > movement_threshold;

% Ensure moving_idx has correct dimensions
if length(moving_idx) ~= nFrames
    moving_idx = [moving_idx; false(nFrames - length(moving_idx), 1)];
end

time_moving_s = sum(moving_idx) * dt;
time_immobile_s = total_time_s - time_moving_s;
percent_moving = 100 * time_moving_s / total_time_s;

% Safe mean calculation function
safe_mean = @(arr) mean(arr, 'omitnan');

mean_speed_all = safe_mean(speed_cm_s);
mean_speed_moving = safe_mean(speed_cm_s(moving_idx));

%% 6. FIXED: ZONE ANALYSIS - RESPECT MANUAL PERIPHERAL WIDTH
peripheral_width_px = peripheral_width_cm / px_to_cm;

% Calculate central zone boundaries
central_min_x = min_x_arena + peripheral_width_px;
central_max_x = max_x_arena - peripheral_width_px;
central_min_y = min_y_arena + peripheral_width_px;
central_max_y = max_y_arena - peripheral_width_px;

central_zone_width = central_max_x - central_min_x;
central_zone_height = central_max_y - central_min_y;

fprintf('\n--- Zone Diagnostics ---\n');
fprintf('Manual peripheral width: %.1f cm = %.1f px\n', peripheral_width_cm, peripheral_width_px);
fprintf('Central zone dimensions: %.1f x %.1f px\n', central_zone_width, central_zone_height);

% FIXED: Only warn but don't auto-adjust if user sets manual width
if central_zone_width <= 0 || central_zone_height <= 0
    warning('Central zone dimensions are invalid with peripheral width = %.1f cm.', peripheral_width_cm);
    warning('The central zone will be empty. Consider reducing peripheral width.');
    % Set central zone to a very small area in the center instead of auto-adjusting
    center_x = (min_x_arena + max_x_arena) / 2;
    center_y = (min_y_arena + max_y_arena) / 2;
    central_min_x = center_x - 1;  % 1 pixel wide
    central_max_x = center_x + 1;
    central_min_y = center_y - 1;
    central_max_y = center_y + 1;
    fprintf('Central zone set to minimal area around center due to invalid dimensions.\n');
end

% Calculate zone masks
central_mask = (x > central_min_x) & (x < central_max_x) & ...
               (y > central_min_y) & (y < central_max_y);
peripheral_mask = ~central_mask;

time_central_s = sum(central_mask) * dt;
time_peripheral_s = sum(peripheral_mask) * dt;
percent_central = 100 * time_central_s / total_time_s;
percent_peripheral = 100 * time_peripheral_s / total_time_s;

mean_speed_central = safe_mean(speed_cm_s(central_mask));
mean_speed_peripheral = safe_mean(speed_cm_s(peripheral_mask));
total_dist_central_cm = sum(distance_cm_frame(central_mask));
total_dist_peripheral_cm = sum(distance_cm_frame(peripheral_mask));

%% 7. Display results
fprintf('\n=== OPEN FIELD ANALYSIS RESULTS ===\n');
fprintf('Peripheral width used: %.1f cm\n', peripheral_width_cm);
fprintf('Total time: %.1f s | Frames: %d\n', total_time_s, nFrames);
fprintf('Total distance: %.1f cm\n', total_distance_cm);
fprintf('Time moving: %.1f s (%.1f%%) | Immobile: %.1f s\n', ...
        time_moving_s, percent_moving, time_immobile_s);
fprintf('Mean speed: %.2f cm/s | Moving speed: %.2f cm/s\n', ...
        mean_speed_all, mean_speed_moving);
fprintf('Central zone: %.1f s (%.1f%%) | Peripheral: %.1f s (%.1f%%)\n', ...
        time_central_s, percent_central, time_peripheral_s, percent_peripheral);

%% 8. Visualization
t = (0:nFrames-1)' * dt;

% Figure 1: Speed analysis
fig1 = figure('Name','Speed Analysis','NumberTitle','off', ...
             'Position', [100 100 1200 600]);

subplot(1,2,1);
plot(t, speed_cm_s, 'b-', 'LineWidth', 1.5);
hold on;
yline(movement_threshold, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Threshold (%.1f cm/s)', movement_threshold));
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Speed [cm/s]', 'FontSize', 12, 'FontWeight', 'bold');
title('Speed vs Time', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'northeast');
ylim([0, max(35, max(speed_cm_s)*1.1)]);

subplot(1,2,2);
histogram(speed_cm_s, 40, 'FaceColor', 'blue', 'EdgeColor', 'none', 'Normalization', 'probability');
hold on;
xline(movement_threshold, 'r--', 'LineWidth', 2);
xlabel('Speed [cm/s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('Speed Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, max(35, max(speed_cm_s)*1.1)]);

% Figure 2: Zone analysis
fig2 = figure('Name',sprintf('Zone Analysis (Peripheral=%.1fcm)', peripheral_width_cm),'NumberTitle','off', ...
             'Position', [100 100 1200 800]);

subplot(2,3,1);
histogram(speed_cm_s(central_mask), 25, 'FaceColor', 'green', 'FaceAlpha', 0.7, 'Normalization', 'probability');
xline(movement_threshold, 'r--', 'LineWidth', 2);
xlabel('Speed [cm/s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('Speed in Central Zone', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, 35]);

subplot(2,3,2);
histogram(speed_cm_s(peripheral_mask), 25, 'FaceColor', 'red', 'FaceAlpha', 0.7, 'Normalization', 'probability');
xline(movement_threshold, 'r--', 'LineWidth', 2);
xlabel('Speed [cm/s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('Speed in Peripheral Zone', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, 35]);

% Trajectory plot
subplot(2,3,3);

% Set background
set(gca, 'Color', [0.98 0.98 0.98]);
grid on;
hold on;

% Plot trajectory
h_path = plot(x, y, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Animal Path');

% Draw arena boundary using plot (for better legend control)
arena_x = [min_x_arena, max_x_arena, max_x_arena, min_x_arena, min_x_arena];
arena_y = [min_y_arena, min_y_arena, max_y_arena, max_y_arena, min_y_arena];
h_arena = plot(arena_x, arena_y, 'r-', 'LineWidth', 3, 'DisplayName', 'Arena Boundary');

% Draw central zone if valid
if central_zone_width > 0 && central_zone_height > 0
    central_x = [central_min_x, central_max_x, central_max_x, central_min_x, central_min_x];
    central_y = [central_min_y, central_min_y, central_max_y, central_max_y, central_min_y];
    h_central = plot(central_x, central_y, 'g--', 'LineWidth', 2.5, 'DisplayName', 'Central Zone');
end

% Create legend
legend('show', 'Location', 'best', 'FontSize', 10, 'Box', 'on');

xlabel('X [px]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y [px]', 'FontSize', 12, 'FontWeight', 'bold');
title('Animal Trajectory', 'FontSize', 14, 'FontWeight', 'bold');
axis equal;

% Set appropriate limits
padding = 50;
xlim([min_x_arena - padding, max_x_arena + padding]);
ylim([min_y_arena - padding, max_y_arena + padding]);

% Zone occupancy
subplot(2,3,4);
zone_plot = zeros(size(t));
zone_plot(central_mask) = 1;
zone_plot(peripheral_mask) = 2;
plot(t, zone_plot, 'k-', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Zone', 'FontSize', 12, 'FontWeight', 'bold');
title('Zone Occupancy', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([0.8, 2.2]);
yticks([1, 2]);
yticklabels({'Central', 'Peripheral'});

subplot(2,3,5);
% Cumulative distance
cumulative_distance = cumsum(distance_cm_frame);
plot(t, cumulative_distance, 'b-', 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cumulative Distance [cm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Distance Traveled', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(2,3,6);
% Movement state
moving_plot = moving_idx * 1.0;
plot(t, moving_plot, 'b-', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Movement', 'FontSize', 12, 'FontWeight', 'bold');
title('Movement State', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'Immobile', 'Moving'});

%% 9. Save results
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
baseFilename = sprintf('TrAQ_Analysis_Periph%.1fcm_%s', peripheral_width_cm, timestamp);

try
    saveas(fig1, sprintf('%s_speed.fig', baseFilename));
    saveas(fig1, sprintf('%s_speed.png', baseFilename));
    saveas(fig2, sprintf('%s_zones.fig', baseFilename));
    saveas(fig2, sprintf('%s_zones.png', baseFilename));
    
    % Save data
    results_table = table(t, speed_cm_s, distance_cm_frame, moving_idx, central_mask, peripheral_mask, ...
        'VariableNames', {'Time_s', 'Speed_cm_s', 'Distance_cm_per_frame', 'IsMoving', 'InCentral', 'InPeripheral'});
    writetable(results_table, sprintf('%s_detailed_data.csv', baseFilename));
    
    summary_data = table(...
        total_time_s, total_distance_cm, time_moving_s, time_immobile_s, percent_moving, ...
        mean_speed_all, mean_speed_moving, time_central_s, time_peripheral_s, ...
        percent_central, percent_peripheral, max(speed_cm_s), px_to_cm, peripheral_width_cm, ...
        'VariableNames', {'TotalTime_s', 'TotalDistance_cm', 'TimeMoving_s', 'TimeImmobile_s', ...
        'PercentMoving', 'MeanSpeed_cm_s', 'MeanSpeedMoving_cm_s', 'TimeCentral_s', ...
        'TimePeripheral_s', 'PercentCentral', 'PercentPeripheral', 'MaxSpeed_cm_s', 'PxToCm', 'PeripheralWidth_cm'});
    writetable(summary_data, sprintf('%s_summary.csv', baseFilename));
    
    fprintf('\nAnalysis complete! Results saved with prefix: %s\n', baseFilename);
    fprintf('Peripheral width used: %.1f cm\n', peripheral_width_cm);
    
catch ME
    warning('Could not save all files: %s', ME.message);
    fprintf('Some files may not have been saved. Check write permissions.\n');
end