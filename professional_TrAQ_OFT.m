function professional_TrAQ_OFT(varargin)
% Robust TrAQ Open Field Analysis: Locomotion + Zones (Peripheral/Central)
%
% This script provides a comprehensive analysis of open field test data,
% including locomotion, zone-based metrics, and detailed visualizations.
% It is designed to be robust, configurable, and easy to use.
%
% Usage:
%   professional_TrAQ_OFT('dataFile', 'Out_Stim2.mat', 'arena_cm', 70, 'peripheral_width_cm', 10)
%
% Parameters:
%   'dataFile' (char) - Name of the .mat file to load (default: 'Out_Stim2.mat')
%   'arena_cm' (numeric) - Real arena dimension in cm (default: 70)
%   'peripheral_width_cm' (numeric) - Width of the peripheral border in cm (default: 10)
%   'movement_threshold' (numeric) - Speed threshold in cm/s to define movement (default: 2.0)
%   'smoothing_window' (numeric) - Size of the moving average window for speed smoothing (default: 5)
%   'histogram_bins' (numeric) - Number of bins for histograms (default: 40)
%
% Example:
%   professional_TrAQ_OFT('dataFile', 'my_tracking_data.mat', 'arena_cm', 50);

    %% --- Configuration ---
    config = parse_inputs(varargin{:});

    %% --- Data Loading ---
    track = load_data(config.dataFile);

    %% --- Data Extraction and Validation ---
    [time, nFrames, fps] = extract_timebase(track);
    [x, y] = extract_positions(track, nFrames);

    %% --- Pixel-to-cm Calibration ---
    [px_to_cm, arena_bounds] = calibrate_pixels(x, y, config.arena_cm);

    %% --- Speed Calculation ---
    [speed_cm_s, distance_cm_frame, total_distance_cm] = calculate_speed(x, y, px_to_cm, 1/fps, config);

    %% --- Movement Analysis ---
    [time_moving_s, time_immobile_s, percent_moving, mean_speed_all, mean_speed_moving, moving_idx] = analyze_movement(speed_cm_s, 1/fps, config);

    %% --- Zone Analysis ---
    [zone_metrics, zone_masks] = analyze_zones(x, y, speed_cm_s, distance_cm_frame, arena_bounds, px_to_cm, config, 1/fps);

    %% --- Display Results ---
    display_results(config, time(end), nFrames, total_distance_cm, time_moving_s, percent_moving, time_immobile_s, mean_speed_all, mean_speed_moving, zone_metrics);

    %% --- Visualization ---
    visualize_results(time, speed_cm_s, x, y, moving_idx, zone_masks, arena_bounds, zone_metrics, config);

    %% --- Save Results ---
    save_results(time, speed_cm_s, distance_cm_frame, moving_idx, zone_masks, zone_metrics, config, total_distance_cm, time_moving_s, time_immobile_s, percent_moving, mean_speed_all, mean_speed_moving);

end

%% --- Helper Functions ---

function config = parse_inputs(varargin)
    % Parses input arguments and sets default values.
    p = inputParser;
    addParameter(p, 'dataFile', 'Out_Stim2.mat', @ischar);
    addParameter(p, 'arena_cm', 70, @isnumeric);
    addParameter(p, 'peripheral_width_cm', 10, @isnumeric);
    addParameter(p, 'movement_threshold', 2.0, @isnumeric);
    addParameter(p, 'smoothing_window', 5, @isnumeric);
    addParameter(p, 'histogram_bins', 40, @isnumeric);
    parse(p, varargin{:});
    config = p.Results;
end

function track = load_data(dataFile)
    % Loads tracking data from the specified file.
    if ~exist(dataFile, 'file')
        error('Data file not found: %s', dataFile);
    end
    loaded_data = load(dataFile);
    if isfield(loaded_data, 'track')
        track = loaded_data.track;
    else
        error('''track'' variable not found in %s.', dataFile);
    end
    fprintf('Loaded %s\n', dataFile);
end

function [time, nFrames, fps] = extract_timebase(track)
    % Extracts time information (frame rate, number of frames) from the track data.
    if isfield(track, 'Time') && ~isempty(track.Time)
        Time = track.Time(:);
        fps = 1 / mean(diff(Time));
        nFrames = numel(Time);
        dt = 1/fps;
        fprintf('Using track.Time: %.2f fps, %d frames\n', fps, nFrames);
    elseif isfield(track, 'FrameRate')
        fps = track.FrameRate;
        dt = 1 / fps;
        if isfield(track, 'Centroid')
            nFrames = size(track.Centroid, 1);
        else
            nFrames = length(track.Time);
        end
    else
        % Fallback if no time information is found
        fps = 30;
        dt = 1/fps;
        if isfield(track, 'Centroid')
            nFrames = size(track.Centroid, 1);
        else
            nFrames = 1000; % Default number of frames
        end
        warning('No time information found. Using default FPS: %d', fps);
    end
    time = (0:nFrames-1)' * dt;
end

function [x, y] = extract_positions(track, nFrames)
    % Extracts and validates position data (x, y coordinates).
    C = [];
    % Check for possible position fields in order of preference
    position_fields = {'Centroid', 'Head', 'Tail', 'Body_Centroid', 'Position'};
    for i = 1:length(position_fields)
        if isfield(track, position_fields{i})
            C = double(track.(position_fields{i}));
            fprintf('Using %s for position data.\n', position_fields{i});
            break;
        end
    end

    if isempty(C)
        error('No suitable position data found.');
    end

    % Remove rows with NaNs
    C(any(isnan(C), 2), :) = [];

    % Handle different coordinate formats (Nx2, 2xN, 3D)
    if size(C, 2) == 2 && size(C, 1) > 2
        % Standard Nx2 format
        x = C(:,1);
        y = C(:,2);
    elseif size(C, 1) == 2 && size(C, 2) > 2
        % Transpose 2xN format
        x = C(1,:)';
        y = C(2,:)';
    elseif ndims(C) == 3
        % Handle 3D arrays from some tracking software
        x = squeeze(C(:,1,1));
        y = squeeze(C(:,1,2));
    else
        error('Unknown coordinate format: %dx%d', size(C,1), size(C,2));
    end

    % Trim position data to match the number of frames
    actual_frames = min(nFrames, length(x));
    x = x(1:actual_frames);
    y = y(1:actual_frames);
    if actual_frames < nFrames
        warning('Position data is shorter than the number of frames. Truncating analysis.');
    end
end

function [px_to_cm, arena_bounds] = calibrate_pixels(x, y, arena_cm)
    % Calibrates pixel coordinates to cm.
    % Add a 10% margin to the arena bounds to avoid edge effects
    x_margin = (max(x) - min(x)) * 0.1;
    y_margin = (max(y) - min(y)) * 0.1;

    min_x_arena = min(x) - x_margin;
    max_x_arena = max(x) + x_margin;
    min_y_arena = min(y) - y_margin;
    max_y_arena = max(y) + y_margin;

    arena_width_px = max_x_arena - min_x_arena;
    arena_height_px = max_y_arena - min_y_arena;

    % Use the larger dimension for calibration to be conservative
    px_to_cm = arena_cm / max(arena_width_px, arena_height_px);

    arena_bounds.min_x = min_x_arena;
    arena_bounds.max_x = max_x_arena;
    arena_bounds.min_y = min_y_arena;
    arena_bounds.max_y = max_y_arena;
end

function [speed_cm_s, distance_cm_frame, total_distance_cm] = calculate_speed(x, y, px_to_cm, dt, config)
    % Calculates speed and distance traveled.
    if length(x) > 1
        dx = diff(x);
        dy = diff(y);
        distance_px_frame = sqrt(dx.^2 + dy.^2);
        distance_cm_frame = [0; distance_px_frame] * px_to_cm;
    else
        distance_cm_frame = zeros(length(x), 1);
    end

    speed_cm_s_raw = distance_cm_frame / dt;

    % Smooth the speed data with a moving average filter
    if length(speed_cm_s_raw) >= config.smoothing_window
        kernel = ones(config.smoothing_window,1)/config.smoothing_window;
        speed_cm_s = conv(speed_cm_s_raw, kernel, 'same');
        % Preserve the original start and end values
        speed_cm_s(1:floor(config.smoothing_window/2)) = speed_cm_s_raw(1:floor(config.smoothing_window/2));
        speed_cm_s(end-floor(config.smoothing_window/2)+1:end) = speed_cm_s_raw(end-floor(config.smoothing_window/2)+1:end);
    else
        speed_cm_s = speed_cm_s_raw;
    end

    total_distance_cm = sum(distance_cm_frame);
end

function [time_moving_s, time_immobile_s, percent_moving, mean_speed_all, mean_speed_moving, moving_idx] = analyze_movement(speed_cm_s, dt, config)
    % Analyzes movement based on a speed threshold.
    moving_idx = speed_cm_s > config.movement_threshold;
    nFrames = length(speed_cm_s);

    if length(moving_idx) ~= nFrames
        moving_idx = [moving_idx; false(nFrames - length(moving_idx), 1)];
    end

    time_moving_s = sum(moving_idx) * dt;
    time_immobile_s = (nFrames * dt) - time_moving_s;
    percent_moving = 100 * time_moving_s / (nFrames * dt);

    safe_mean = @(arr) mean(arr, 'omitnan');

    mean_speed_all = safe_mean(speed_cm_s);
    mean_speed_moving = safe_mean(speed_cm_s(moving_idx));
end

function [zone_metrics, zone_masks] = analyze_zones(x, y, speed_cm_s, distance_cm_frame, arena_bounds, px_to_cm, config, dt)
    % Performs zone analysis (central vs. peripheral).
    peripheral_width_px = config.peripheral_width_cm / px_to_cm;

    central_min_x = arena_bounds.min_x + peripheral_width_px;
    central_max_x = arena_bounds.max_x - peripheral_width_px;
    central_min_y = arena_bounds.min_y + peripheral_width_px;
    central_max_y = arena_bounds.max_y - peripheral_width_px;

    central_zone_width = central_max_x - central_min_x;
    central_zone_height = central_max_y - central_min_y;

    % Handle cases where the peripheral width is too large
    if central_zone_width <= 0 || central_zone_height <= 0
        warning('Central zone dimensions are invalid. The central zone will be empty.');
        center_x = (arena_bounds.min_x + arena_bounds.max_x) / 2;
        center_y = (arena_bounds.min_y + arena_bounds.max_y) / 2;
        central_min_x = center_x - 1;
        central_max_x = center_x + 1;
        central_min_y = center_y - 1;
        central_max_y = center_y + 1;
    end

    central_mask = (x > central_min_x) & (x < central_max_x) & ...
                   (y > central_min_y) & (y < central_max_y);
    peripheral_mask = ~central_mask;

    zone_masks.central = central_mask;
    zone_masks.peripheral = peripheral_mask;

    total_time_s = length(x) * dt;
    zone_metrics.time_central_s = sum(central_mask) * dt;
    zone_metrics.time_peripheral_s = sum(peripheral_mask) * dt;
    zone_metrics.percent_central = 100 * zone_metrics.time_central_s / total_time_s;
    zone_metrics.percent_peripheral = 100 * zone_metrics.time_peripheral_s / total_time_s;

    safe_mean = @(arr) mean(arr, 'omitnan');
    zone_metrics.mean_speed_central = safe_mean(speed_cm_s(central_mask));
    zone_metrics.mean_speed_peripheral = safe_mean(speed_cm_s(peripheral_mask));
    zone_metrics.total_dist_central_cm = sum(distance_cm_frame(central_mask));
    zone_metrics.total_dist_peripheral_cm = sum(distance_cm_frame(peripheral_mask));

    zone_metrics.central_bounds.min_x = central_min_x;
    zone_metrics.central_bounds.max_x = central_max_x;
    zone_metrics.central_bounds.min_y = central_min_y;
    zone_metrics.central_bounds.max_y = central_max_y;
end

function display_results(config, total_time_s, nFrames, total_distance_cm, time_moving_s, percent_moving, time_immobile_s, mean_speed_all, mean_speed_moving, zone_metrics)
    % Displays analysis results in the command window.
    fprintf('\n=== OPEN FIELD ANALYSIS RESULTS ===\n');
    fprintf('Peripheral width used: %.1f cm\n', config.peripheral_width_cm);
    fprintf('Total time: %.1f s | Frames: %d\n', total_time_s, nFrames);
    fprintf('Total distance: %.1f cm\n', total_distance_cm);
    fprintf('Time moving: %.1f s (%.1f%%) | Immobile: %.1f s\n', ...
            time_moving_s, percent_moving, time_immobile_s);
    fprintf('Mean speed: %.2f cm/s | Moving speed: %.2f cm/s\n', ...
            mean_speed_all, mean_speed_moving);
    fprintf('Central zone: %.1f s (%.1f%%) | Peripheral: %.1f s (%.1f%%)\n', ...
            zone_metrics.time_central_s, zone_metrics.percent_central, zone_metrics.time_peripheral_s, zone_metrics.percent_peripheral);
end

function visualize_results(t, speed_cm_s, x, y, moving_idx, zone_masks, arena_bounds, zone_metrics, config)
    % Creates and saves figures for visualization.

    % --- Figure 1: Speed Analysis ---
    fig1 = figure('Name','Speed Analysis','NumberTitle','off', ...
                 'Position', [100 100 1200 600]);

    ax1 = subplot(1,2,1);
    plot(t, speed_cm_s, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(config.movement_threshold, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Threshold (%.1f cm/s)', config.movement_threshold));
    xlabel('Time [s]');
    ylabel('Speed [cm/s]');
    title('Speed vs Time');
    grid on;
    legend('show', 'Location', 'northeast');
    ylim([0, max(speed_cm_s)*1.1]);
    set(ax1, 'FontSize', 12, 'FontWeight', 'bold');

    ax2 = subplot(1,2,2);
    histogram(speed_cm_s, config.histogram_bins, 'FaceColor', 'blue', 'EdgeColor', 'none', 'Normalization', 'probability');
    hold on;
    xline(config.movement_threshold, 'r--', 'LineWidth', 2);
    xlabel('Speed [cm/s]');
    ylabel('Probability');
    title('Speed Distribution');
    grid on;
    xlim([0, max(speed_cm_s)*1.1]);
    set(ax2, 'FontSize', 12, 'FontWeight', 'bold');

    % --- Figure 2: Zone Analysis ---
    fig2 = figure('Name',sprintf('Zone Analysis (Peripheral=%.1fcm)', config.peripheral_width_cm),'NumberTitle','off', ...
                 'Position', [100 100 1200 800]);

    % Speed in central zone
    ax3 = subplot(2,3,1);
    histogram(speed_cm_s(zone_masks.central), config.histogram_bins, 'FaceColor', 'green', 'FaceAlpha', 0.7, 'Normalization', 'probability');
    xline(config.movement_threshold, 'r--', 'LineWidth', 2);
    xlabel('Speed [cm/s]');
    ylabel('Probability');
    title('Speed in Central Zone');
    grid on;
    xlim([0, max(speed_cm_s)*1.1]);
    set(ax3, 'FontSize', 12, 'FontWeight', 'bold');

    % Speed in peripheral zone
    ax4 = subplot(2,3,2);
    histogram(speed_cm_s(zone_masks.peripheral), config.histogram_bins, 'FaceColor', 'red', 'FaceAlpha', 0.7, 'Normalization', 'probability');
    xline(config.movement_threshold, 'r--', 'LineWidth', 2);
    xlabel('Speed [cm/s]');
    ylabel('Probability');
    title('Speed in Peripheral Zone');
    grid on;
    xlim([0, max(speed_cm_s)*1.1]);
    set(ax4, 'FontSize', 12, 'FontWeight', 'bold');

    % Trajectory plot
    ax5 = subplot(2,3,3);
    set(ax5, 'Color', [0.98 0.98 0.98]);
    grid on;
    hold on;
    plot(x, y, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Animal Path');
    arena_x = [arena_bounds.min_x, arena_bounds.max_x, arena_bounds.max_x, arena_bounds.min_x, arena_bounds.min_x];
    arena_y = [arena_bounds.min_y, arena_bounds.min_y, arena_bounds.max_y, arena_bounds.max_y, arena_bounds.min_y];
    plot(arena_x, arena_y, 'r-', 'LineWidth', 3, 'DisplayName', 'Arena Boundary');

    central_bounds = zone_metrics.central_bounds;
    if (central_bounds.max_x - central_bounds.min_x) > 0 && (central_bounds.max_y - central_bounds.min_y) > 0
        central_x = [central_bounds.min_x, central_bounds.max_x, central_bounds.max_x, central_bounds.min_x, central_bounds.min_x];
        central_y = [central_bounds.min_y, central_bounds.min_y, central_bounds.max_y, central_bounds.max_y, central_bounds.min_y];
        plot(central_x, central_y, 'g--', 'LineWidth', 2.5, 'DisplayName', 'Central Zone');
    end

    legend('show', 'Location', 'best', 'FontSize', 10, 'Box', 'on');
    xlabel('X [px]');
    ylabel('Y [px]');
    title('Animal Trajectory');
    axis equal;
    padding = 50;
    xlim([arena_bounds.min_x - padding, arena_bounds.max_x + padding]);
    ylim([arena_bounds.min_y - padding, arena_bounds.max_y + padding]);
    set(ax5, 'FontSize', 12, 'FontWeight', 'bold');

    % Zone occupancy
    ax6 = subplot(2,3,4);
    zone_plot = zeros(size(t));
    zone_plot(zone_masks.central) = 1;
    zone_plot(zone_masks.peripheral) = 2;
    plot(t, zone_plot, 'k-', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Zone');
    title('Zone Occupancy');
    grid on;
    ylim([0.8, 2.2]);
    yticks([1, 2]);
    yticklabels({'Central', 'Peripheral'});
    set(ax6, 'FontSize', 12, 'FontWeight', 'bold');

    % Cumulative distance
    ax7 = subplot(2,3,5);
    cumulative_distance = cumsum(distance_cm_frame);
    plot(t, cumulative_distance, 'b-', 'LineWidth', 2);
    xlabel('Time [s]');
    ylabel('Cumulative Distance [cm]');
    title('Distance Traveled');
    grid on;
    set(ax7, 'FontSize', 12, 'FontWeight', 'bold');

    % Movement state
    ax8 = subplot(2,3,6);
    moving_plot = moving_idx * 1.0;
    plot(t, moving_plot, 'b-', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Movement');
    title('Movement State');
    grid on;
    ylim([-0.1, 1.1]);
    yticks([0, 1]);
    yticklabels({'Immobile', 'Moving'});
    set(ax8, 'FontSize', 12, 'FontWeight', 'bold');

    timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
    baseFilename = sprintf('TrAQ_Analysis_Periph%.1fcm_%s', config.peripheral_width_cm, timestamp);
    try
        saveas(fig1, sprintf('%s_speed.png', baseFilename));
        saveas(fig2, sprintf('%s_zones.png', baseFilename));
    catch ME
        warning('Could not save figures: %s', ME.message);
    end
end

function save_results(t, speed_cm_s, distance_cm_frame, moving_idx, zone_masks, zone_metrics, config, total_distance_cm, time_moving_s, time_immobile_s, percent_moving, mean_speed_all, mean_speed_moving)
    % Saves analysis results to CSV files.
    timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
    baseFilename = sprintf('TrAQ_Analysis_Periph%.1fcm_%s', config.peripheral_width_cm, timestamp);

    try
        results_table = table(t, speed_cm_s, distance_cm_frame, moving_idx, zone_masks.central, zone_masks.peripheral, ...
            'VariableNames', {'Time_s', 'Speed_cm_s', 'Distance_cm_per_frame', 'IsMoving', 'InCentral', 'InPeripheral'});
        writetable(results_table, sprintf('%s_detailed_data.csv', baseFilename));

        summary_data = table(...
            t(end), total_distance_cm, time_moving_s, time_immobile_s, percent_moving, ...
            mean_speed_all, mean_speed_moving, zone_metrics.time_central_s, zone_metrics.time_peripheral_s, ...
            zone_metrics.percent_central, zone_metrics.percent_peripheral, max(speed_cm_s), config.peripheral_width_cm, ...
            'VariableNames', {'TotalTime_s', 'TotalDistance_cm', 'TimeMoving_s', 'TimeImmobile_s', ...
            'PercentMoving', 'MeanSpeed_cm_s', 'MeanSpeedMoving_cm_s', 'TimeCentral_s', ...
            'TimePeripheral_s', 'PercentCentral', 'PercentPeripheral', 'MaxSpeed_cm_s', 'PeripheralWidth_cm'});
        writetable(summary_data, sprintf('%s_summary.csv', baseFilename));

        fprintf('\nAnalysis complete! Results saved with prefix: %s\n', baseFilename);
    catch ME
        warning('Could not save all files: %s', ME.message);
    end
end
