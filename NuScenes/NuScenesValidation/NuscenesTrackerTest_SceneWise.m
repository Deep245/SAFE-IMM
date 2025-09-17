
clc; clear;
run("varInit.m"); % expects window_size, clusterer, thresholds, etc.

%% === Load scene .mat file ===
scene   = load('nuscenes_scene_mat_output/scene-0103.mat');
samples = scene.samples;
num_frames = length(samples);

% === Assign consistent color per object (track) ===
all_ids = {};
for i = 1:num_frames
    for j = 1:length(samples{i}.gt_boxes)
        all_ids{end+1} = samples{i}.gt_boxes{j}.instance_token; %#ok<AGROW>
    end
end
unique_ids = unique(all_ids);
colors     = hsv(length(unique_ids));
id_map     = containers.Map(unique_ids, mat2cell(colors, ones(1,length(unique_ids)), 3));

% === Track histories and clustering state ===
tracks = containers.Map;
cluster_history_static  = cell(1, num_frames);
cluster_history_dynamic = cell(1, num_frames);
idxS = []; idxD = [];
static_points_window  = [];
dynamic_points_window = [];


% === Time index (seconds; name kept for API compatibility) ===
unique_microseconds = [];  % actually seconds; kept to match clustering function signature


% === Figure ===
figure('Name', 'nuScenes Top-down Visualization (3 Views)', 'NumberTitle', 'off');

trackerTime = 0;                          % numeric scalar
prev_tstamp = double(samples{1}.timestamp) * 1e-6;
if ~exist('trackHist','var') || isempty(trackHist)
    trackHist = containers.Map('KeyType','double','ValueType','any'); % id -> [x y; ...]
end
histLen = 50;   % how many past points to keep per track

% ===== SELECT IMM VARIANT TO TEST =====
testVariant = "linear";   % "alpabetagamma" | "linear" | "nonlinear" | "baseline" 

switch testVariant
    case "linear"
        initFcn = @initIMM_CV_CA_CT_linear;
    case "nonlinear"
        initFcn = @initIMM_CV_CA_CT_nonlinear;
    case "baseline"
        initFcn = @initIMM_CV_CA_CT_baseline;
    otherwise
        error('Unknown testVariant: %s', testVariant);
end
posStd = 1.5;  velStd = 11.0;       % measurement std devs, this tuned for newIMM
%posStd = 1;  velStd = 4.0;       % measurement std devs, this is tuned for
% Conventional IMM
R = diag([posStd^2 posStd^2 posStd^2 velStd^2 velStd^2 velStd^2]);
coastFrames  = 3;    % allow prediction-only up to 8 frames
% --- Build the trackerGNN that uses the IMM initializer ---
trk = trackerGNN( ...
    'FilterInitializationFcn', initFcn, ...
    'AssignmentThreshold', 30, ...
    'ConfirmationThreshold', [2 5], ...
    'DeletionThreshold', coastFrames);
timePrevCost  = [];     % the time those 


% % ===== PROFILING BUFFERS =====
% t_frame_total   = zeros(num_frames,1);  % wall time per frame (s) incl. plotting
% t_tracker_only  = zeros(num_frames,1);  % time inside trackerGNN (s)
% fps_inst        = zeros(num_frames,1);  % instantaneous FPS = 1/t_frame_total
% fps_rolling     = zeros(num_frames,1);  % rolling FPS over last ROLL frames
% ndets_per_frame = zeros(num_frames,1);  % detections fed to tracker
% ntracks_out     = zeros(num_frames,1);  % confirmed tracks returned
% state_dim_hist  = zeros(num_frames,1);  % median track state length (6 CV, 9 CA, 7 CT, 0 none)
% ROLL = 30;                               % rolling window length
% doPlot = true;                           % set false for pure algorithm speed

for i = 1:num_frames
    s = samples{i};

    % === Time (seconds, consistent) ===
    t = double(s.timestamp) * 1e-6;
    unique_microseconds = [unique_microseconds; t]; %#ok<AGROW>

    tstamp = double(s.timestamp) * 1e-6;    
    if i == 1
        trackerTime = 0;                      % stay numeric
    else
        Ts = max(tstamp - prev_tstamp, 1e-3); % positive step
        trackerTime = trackerTime + Ts;       % stays double
    end
    prev_tstamp = tstamp;    
    % defensive: if it somehow became a cell, unwrap
    if iscell(trackerTime), trackerTime = trackerTime{1}; end
    trackerTime = double(trackerTime);


    % === Transforms (compute once) ===
    R_radar = quat2rotm(s.radar_rotation);   % [w x y z]
    t_radar = s.radar_translation(:);

    R_ego   = quat2rotm(s.ego_rotation);
    t_ego   = s.ego_translation(:);

    R_global = R_ego * R_radar;
    t_global = R_ego * t_radar + t_ego;

    % Global yaw of radar x-axis (for oriented 2D boxes)
    yaw_global = atan2(R_global(2,1), R_global(1,1));

    % === Radar points (positions) ===
    radar_pts_radar = s.radar_points(:, 1:3)';   % 3xN
    pts_ego    = R_radar * radar_pts_radar + t_radar;
    pts_global = R_ego   * pts_ego + t_ego;

    % === Radar velocities (rotate only; do NOT translate) ===
    v_radar  = s.radar_points(:,4:6)';   % 3xN
    v_ego    = R_radar * v_radar;
    v_global = R_ego   * v_ego;


    % === Predefine outputs for safety ===
    staticDetection = []; static_yaws = [];
    static_global   = []; dynamic_global = [];

    % === Bounds & detection table (x y z vx vy vz) ===
    radarPoints = s.radar_points(:, 1:6); % x,y,z,vx,vy,vz (radar frame)
    xLowerBound = 0;   xHigerBound = 80;
    yLowerBound = -20; yHigerBound = 20;
    zLowerBound = -1;  zHigerBound = 1;
    detections = processDet(radarPoints, t, ...
        xLowerBound, xHigerBound, yLowerBound, yHigerBound, zLowerBound, zHigerBound);

    % === Build objectDetection array (Time in seconds) ===
    ObjDet = cell(1, height(detections));
    for l = 1:height(detections)
        ObjDet{l} = objectDetection( ...
            detections.t(l), ...
            [detections.X(l), detections.Y(l), detections.Z(l), ...
             detections.Vx(l), detections.Vy(l), detections.Vz(l)], ...
            "MeasurementNoise", measurement_noise, "MeasurementParameters", mp);
    end

    % === Static / Dynamic classification ===
    [staticmeas, dynamicmeas, statictime, dynamictime] = ClassifyDet(ObjDet);

    % === Transform classified points to GLOBAL for plotting ===
    if ~isempty(staticmeas)
        static_pts_radar = staticmeas(:, 1:3)';          % 3xNs
        spts_ego    = R_radar * static_pts_radar + t_radar;
        static_global = R_ego * spts_ego + t_ego;        % 3xNs
    end
    if ~isempty(dynamicmeas)
        dynamic_pts_radar = dynamicmeas(:, 1:3)';        % 3xNd
        dpts_ego    = R_radar * dynamic_pts_radar + t_radar;
        dynamic_global = R_ego * dpts_ego + t_ego;       % 3xNd
    end

    % === Clustering (sliding window over time indices) ===
    [idxD, idxS, cluster_history_static, cluster_history_dynamic] = ...
        clusteringRadarPoints(i, window_size, statictime, dynamictime, ...
                              staticmeas, dynamicmeas, unique_microseconds, ...
                              clusterer, idxD, idxS, ...
                              dynamic_points_window, static_points_window, ...
                              cluster_history_static, cluster_history_dynamic);


    % % Static detection
    % if ~isempty(cluster_history_static{i})
    %     [cluster_history_static,staticDetection] = update_staticCandiate(i,cluster_history_static,first_satisfied_index,persistence_threshold,threshold_static,...
    %     merge_threshold_x,merge_threshold_y,merge_threshold_z);
    %     boxes_g = transform_boxes_to_global(staticDetection, R_global, t_global);
    %     static_yaws = repmat(yaw_global, 1, size(boxes_g, 1));
    % end

    % Dynamic detection (writes back D,L,pbb via assignin inside the function)
    if ~isempty(cluster_history_dynamic{i})
        [cluster_history_dynamic,DynamicDetection,pbb_c] = update_dynamicCandidates(i,cluster_history_dynamic,dynamic_threshold,DynamicBoxThreshold,...
        merge_threshold_x,merge_threshold_y,merge_threshold_z,frames_dBox);
        DynamicDetection = transform_boxes_to_global(DynamicDetection, R_global, t_global);
        dynamic_yaws = repmat(yaw_global, 1, size(DynamicDetection, 1));
        
        pbb_global = to_global_pbb(pbb_c, R_global, t_global);
        % pbb_global = [pbb_global(:,1:3),pbb_c(:,4:6)];
        v_global_pbb = (R_global * pbb_c(:,4:6).').';  % rotate velocities to GLOBAL (no translate)
        pbb_global = [pbb_global(:,1:3), v_global_pbb]; % pos+vel now in the SAME frame
        
        dets = cell(1, size(pbb_global,1));
        for j = 1:size(pbb_global,1)
            meas = pbb_global(j,:).';                   % 6x1 vector
            dets{j} = objectDetection(trackerTime, meas, 'MeasurementNoise', R);
        end

        % — also expose IDs for plotting (use current labels) —
        cluster_history_dynamic{i}.ids_for_plot = cluster_history_dynamic{i}.dynamic_box_label(:);
    end
    
    % trac = trk(dets, trackerTime);
    % ===== time tracker =====
    ndets = 0;
    if exist('dets','var') && ~isempty(dets), ndets = numel(dets); end
    ndets_per_frame(i) = ndets;    
    tic_total = tic;    
    tic_trk = tic;
    trac = trk(dets, trackerTime);
    t_tracker_only(i) = toc(tic_trk);    
    % keep only confirmed & aged (your existing post-filter)
    if ~isempty(trac)
        trac = trac([trac.IsConfirmed]);
        ages = [trac.Age];
        trac = trac(ages > 3);
    end
    ntracks_out(i) = numel(trac);    
    % state-dimension proxy (which model dominates)
    if ~isempty(trac)
        stlens = arrayfun(@(tt) numel(tt.State), trac);
        state_dim_hist(i) = median(stlens);
    else
        state_dim_hist(i) = 0;
    end    
    % % ===== plotting toggle =====
    % if ~doPlot
    %     % skip the plotting block entirely by jumping to end-of-loop label
    %     t_frame_total(i) = toc(tic_total);
    %     fps_inst(i)      = 1 / max(t_frame_total(i), eps);
    %     lo = max(1, i-ROLL+1);
    %     fps_rolling(i)   = (i-lo+1) / sum(t_frame_total(lo:i));
    %     continue;   % <-- bypass plots if disabled
    % end

    
    % keep only confirmed for plotting (optional but recommended)
    if ~isempty(trac)
        trac = trac([trac.IsConfirmed]);
        ages = [trac.Age];
        trac = trac(ages > 3);
    end
    

    % % ===== finalize per-frame timing (if plotting enabled) =====
    % t_frame_total(i) = toc(tic_total);
    % fps_inst(i)      = 1 / max(t_frame_total(i), eps);
    % lo = max(1, i-ROLL+1);
    % fps_rolling(i)   = (i-lo+1) / sum(t_frame_total(lo:i));



    %% === Subplot 1: Full Map + Radar + GT Boxes + Tracks ===
    subplot(1,3,1); cla; hold on; axis equal; grid on;
    title(sprintf('Top-down Map View | Sample %d / %d', i, num_frames));
    xlabel('Global X (m)'); ylabel('Global Y (m)');

    % Plot map lanes
    for j = 1:length(s.lanes)
        lane_poly = s.lanes{j};
        if size(lane_poly,2) == 2
            plot(lane_poly(:,1), lane_poly(:,2), 'k-', 'LineWidth', 1);
        end
    end

    % Plot radar points (GLOBAL)
    scatter(pts_global(1,:), pts_global(2,:), 8, 'r', 'filled');

    % Ego pose + heading in GLOBAL
    plot(t_ego(1), t_ego(2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    heading_len = 3;
    dir_vec = R_ego * [heading_len; 0; 0];
    quiver(t_ego(1), t_ego(2), dir_vec(1), dir_vec(2), 0, 'b', 'LineWidth', 2);


    % ---- front-sector gates ----
    maxRange    = 100;      % meters
    fovRadarDeg = 60;       % total FoV for radar front (±30°)
    fovCamDeg   = 60;       % total FoV for camera front (±35°)

    % radar forward yaw in GLOBAL (radar x-axis after ego transform)
    yawRadar = atan2(R_global(2,1), R_global(1,1));

    % ego forward yaw in GLOBAL (ego x-axis)
    yawEgo = atan2(R_ego(2,1), R_ego(1,1));

    % camera yaw: try to read it; otherwise assume camera aligned with ego forward
    if isfield(s, 'camera_rotation')
        R_cam     = quat2rotm(s.camera_rotation);   % expected [w x y z]
        R_cam_glo = R_ego * R_cam;
        yawCam    = atan2(R_cam_glo(2,3), R_cam_glo(1,3));  % camera z-axis forward
    else
        % FALLBACK: no camera extrinsics → use ego heading as camera forward
        yawCam = yawEgo;
    end

    halfFovRadar = deg2rad(fovRadarDeg/2);
    halfFovCam   = deg2rad(fovCamDeg/2);



    % Ground-truth boxes + track trails
    % Ground-truth boxes + track trails (front radar/camera only)
    for j = 1:length(s.gt_boxes)
        box   = s.gt_boxes{j};
        center = box.center(:);             % [x;y;z] GLOBAL
        size3  = box.size(:);               % [w,l,h]
        quat   = box.rotation(:);
        R_box  = quat2rotm(quat');

        % --- gating in XY (range & angle) ---
        dx = center(1) - t_ego(1);
        dy = center(2) - t_ego(2);
        rngXY = hypot(dx, dy);
        if rngXY > maxRange, continue; end

        theta    = atan2(dy, dx);
        dYawRad  = atan2(sin(theta - yawRadar), cos(theta - yawRadar));
        dYawCam  = atan2(sin(theta - yawCam),   cos(theta - yawCam));

        inRadar  = abs(dYawRad) <= halfFovRadar;
        inCamera = abs(dYawCam) <= halfFovCam;

        % pick one of these:
        showIt = (inRadar || inCamera);   % union (front radar OR front camera)
        % showIt = inRadar;               % radar-front only (use this if no camera)
        % showIt = (inRadar && inCamera); % intersection (only if in both)

        if ~showIt, continue; end


        % --- draw oriented 2D footprint on XY plane ---
        l = size3(2);  % length along box local x
        w = size3(1);  % width  along box local y
        local_corners = 0.5 * [l w; -l w; -l -w; l -w]';
        global_corners = R_box(1:2,1:2) * local_corners + center(1:2);

        % color by instance id; edge color shows which sensor(s) gated it in
        color = id_map(box.instance_token);
        % NEW
        if inRadar && inCamera
            edgeCol = [0 0 0];          % both
        elseif inRadar
            edgeCol = [0 0.6 0];        % radar-only
        elseif inCamera
            edgeCol = [0 0 0.8];        % camera-only
        else
            continue
        end

        fill(global_corners(1,:), global_corners(2,:), color, ...
             'FaceAlpha', 0.5, 'EdgeColor', edgeCol, 'LineWidth', 1.5);

        % optional: only keep history for shown boxes
        if isKey(tracks, box.instance_token)
            track = tracks(box.instance_token);
            track(end+1,:) = center(1:2)'; %#ok<AGROW>
            tracks(box.instance_token) = track;
        else
            tracks(box.instance_token) = center(1:2)';
        end
    end


    keys_ = tracks.keys;
    for k = 1:length(keys_)
        track = tracks(keys_{k});
        plot(track(:,1), track(:,2), '-', 'Color', id_map(keys_{k}), 'LineWidth', 1);
    end

    xlim([t_ego(1) - 40, t_ego(1) + 40]);
    ylim([t_ego(2) - 40, t_ego(2) + 40]);

    %% === Subplot 2: Camera View ===
    subplot(1,3,2); cla;
    imshow(uint8(s.camera_image_data));
    title('Camera Front');

    %% === Subplot 3: Clean Map + Radar + Ego + Static/Dynamic ===
    subplot(1,3,3); cla; hold on; axis equal; grid on;
    title('Radar + Ego Only');
    xlabel('Global X (m)'); ylabel('Global Y (m)');

    % Plot only map lanes
    for j = 1:length(s.lanes)
        lane_poly = s.lanes{j};
        if size(lane_poly,2) == 2
            plot(lane_poly(:,1), lane_poly(:,2), 'k-', 'LineWidth', 1);
        end
    end

    % Static (blue) and dynamic (red) points in GLOBAL
    if ~isempty(static_global)
        scatter(static_global(1,:), static_global(2,:), 20, 'b', 'filled');
    end
    if ~isempty(dynamic_global)
        scatter(dynamic_global(1,:), dynamic_global(2,:), 20, 'r', 'filled');
    end

    % % Oriented 2D bounding boxes for persistent static clusters
    % if ~isempty(boxes_g)
    %     plot_bounding_boxes_2d(boxes_g, static_yaws, 'g', 2);
    % end
    if ~isempty(DynamicDetection)
        plot_bounding_boxes_2d(DynamicDetection, dynamic_yaws, 'r', 2);
    end

    % --- Update history with current tracks ---
    if exist('trac','var') && ~isempty(trac)
        idsNow = [trac.TrackID];

        for it = 1:numel(trac)
            id = trac(it).TrackID;

            [pRow, vRow] = posvelFromTrackCompat(trac(it)); % <-- new
            px = pRow(1); py = pRow(2); vx = vRow(1); vy = vRow(2);
            

            % append to history
            if ~exist('trackHist','var') || isempty(trackHist)
                trackHist = containers.Map('KeyType','double','ValueType','any');
            end
            if isKey(trackHist, id), H = trackHist(id); else, H = zeros(0,2); end

            H = [H; [px py]];
            histLen = 50;
            if size(H,1) > histLen, H = H(end-histLen+1:end,:); end
            trackHist(id) = H;
        end

        % draw histories
        klist = trackHist.keys;
        NEAR_ZERO = 0.2;     % tune to your map scale
        MAX_JUMP  = 20;      % tune to your scene dynamics
        
        for k = 1:numel(klist)
            idk = klist{k};
            H = trackHist(idk);
            if isempty(H), continue; end
            Hs = sanitizeHistoryForPlot(H, NEAR_ZERO, MAX_JUMP);
            plot(Hs(:,1), Hs(:,2), 'g-', 'LineWidth', 1.5);
        end

        % draw current states
        velScale = 0.3;  % tune arrow length
        for it = 1:numel(trac)
            id = trac(it).TrackID;
            [pRow, vRow] = posvelFromTrackCompat(trac(it));
            plot(pRow(1), pRow(2), 'go', 'MarkerSize', 6, 'LineWidth', 1.5);
            quiver(pRow(1), pRow(2), velScale*vRow(1), velScale*vRow(2), 0, 'g', 'LineWidth', 1);
            text(pRow(1)+0.5, pRow(2)+0.5, sprintf('%d', id), 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
        end

        % (optional) drop history for deleted tracks
        gone = setdiff(cell2mat(trackHist.keys), idsNow);
        %for g = gone, remove(trackHist, g); end
    end



    % Ego pose + heading
    plot(t_ego(1), t_ego(2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    quiver(t_ego(1), t_ego(2), dir_vec(1), dir_vec(2), 0, 'b', 'LineWidth', 2);

    % Optional: velocity arrows of raw radar points (GLOBAL)
    % quiver(pts_global(1,:), pts_global(2,:), vx_global, vy_global, 0, 'm', 'LineWidth', 1);

    xlim([t_ego(1) - 40, t_ego(1) + 40]);
    ylim([t_ego(2) - 40, t_ego(2) + 40]);

    drawnow limitrate;
end


% % ===== SUMMARY & SAVE =====
% valid = t_frame_total > 0;
% mean_fps   = mean(fps_inst(valid));
% p50_time   = prctile(t_frame_total(valid), 50);
% p95_time   = prctile(t_frame_total(valid), 95);
% mean_trk   = mean(t_tracker_only(valid));
% p95_trk    = prctile(t_tracker_only(valid), 95);
% 
% fprintf('\n=== Performance [%s] ===\n', testVariant);
% fprintf('Frames: %d\n', nnz(valid));
% fprintf('Mean FPS (total frame): %.1f\n', mean_fps);
% fprintf('Median frame time: %.3f s | 95th: %.3f s\n', p50_time, p95_time);
% fprintf('Tracker time: mean %.3f s | 95th %.3f s\n', mean_trk, p95_trk);
% fprintf('Mean detections/frame: %.1f | Mean confirmed tracks: %.1f\n', ...
%     mean(ndets_per_frame(valid)), mean(ntracks_out(valid)));
% fprintf('State length median counts 6:%d 7:%d 9:%d 0:%d\n', ...
%     nnz(state_dim_hist==6), nnz(state_dim_hist==7), nnz(state_dim_hist==9), nnz(state_dim_hist==0));
% 
% T = table( (1:num_frames).', ndets_per_frame, t_tracker_only, t_frame_total, ...
%            fps_inst, fps_rolling, state_dim_hist, ...
%            'VariableNames', {'frame','ndets','t_tracker','t_total','fps_inst','fps_roll','state_dim'});
% csvname = sprintf('perf_%s_%s.csv', testVariant, datestr(now,'yyyymmdd_HHMMSS'));
% writetable(T, csvname);
% fprintf('Saved per-frame metrics to %s\n', csvname);




% === Helper: draw oriented 2D boxes ===
function plot_bounding_boxes_2d(bounding_boxes, yaws, color, linewidth)
    if nargin < 3, color = 'r'; end
    if nargin < 4, linewidth = 2; end
    for i = 1:size(bounding_boxes, 1)
        box = bounding_boxes(i, :);
        if numel(box) ~= 6, continue; end
        xmin = box(1); ymin = box(2); xmax = box(4); ymax = box(5);
        width  = xmax - xmin;  height = ymax - ymin;
        if width <= 0 || height <= 0, continue; end
        cx = (xmin + xmax) / 2;  cy = (ymin + ymax) / 2;
        angle_rad = (i <= numel(yaws)) * yaws(i);  % 0 if missing
        rect_x = [-width/2, width/2, width/2, -width/2, -width/2];
        rect_y = [-height/2, -height/2, height/2,  height/2, -height/2];
        R = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
        rotated = R * [rect_x; rect_y];
        plot(rotated(1,:) + cx, rotated(2,:) + cy, color, 'LineWidth', linewidth);
    end
end

% function imm = initIMM_CV_CA_CT_linear(detection)
% % IMM with CV (KF), CA (KF), and CT (EKF)
% % detection.Measurement = [x y z vx vy vz], detection.MeasurementNoise = 6x6
% 
% R6 = detection.MeasurementNoise;
% m  = detection.Measurement(:);   % force column
% 
% x0  = m(1);  y0  = m(2);  z0  = m(3);
% vx0 = (numel(m)>=4)*m(4);
% vy0 = (numel(m)>=5)*m(5);
% vz0 = (numel(m)>=6)*m(6);
% 
% %% ===== CV (trackingKF), state: [x; vx; y; vy; z; vz] (6x1)
% kfCV = trackingKF( ...
%     'MotionModel','3D Constant Velocity', ...
%     'State',[x0; vx0; y0; vy0; z0; vz0], ...      % <-- column
%     'StateCovariance',diag([4 25 4 25 4 25]), ...
%     'MeasurementModel',eye(6), ...
%     'MeasurementNoise',R6);
% 
% %% ===== CA (trackingKF), state: [x; vx; ax; y; vy; ay; z; vz; az] (9x1)
% Hca = [ ...
%     1 0 0   0 0 0   0 0 0;   % x
%     0 0 0   1 0 0   0 0 0;   % y
%     0 0 0   0 0 0   1 0 0;   % z
%     0 1 0   0 0 0   0 0 0;   % vx
%     0 0 0   0 1 0   0 0 0;   % vy
%     0 0 0   0 0 0   0 1 0];  % vz
% 
% kfCA = trackingKF( ...
%     'MotionModel','3D Constant Acceleration', ...
%     'State',[x0; vx0; 0;  y0; vy0; 0;  z0; vz0; 0], ...   % <-- column
%     'StateCovariance',diag([4 25 25  4 25 25  9 36 36]), ...
%     'MeasurementModel',Hca, ...
%     'MeasurementNoise',R6);
% 
% %% ===== CT (trackingEKF), state: [x; vx; y; vy; omega; z; vz] (7x1)
% xct0 = [x0; vx0; y0; vy0; 0; z0; vz0];
% Pct0 = diag([4 25 4 25 (pi/2)^2 9 36]);
% Qct  = blkdiag(0.2*eye(2),0.2*eye(2),0.02,0.4*eye(2));
% measFcnCT = @(x)[x(1); x(3); x(6); x(2); x(4); x(7)];  % [x y z vx vy vz]
% 
% kfCT = trackingEKF( ...
%     'State',xct0, ...
%     'StateCovariance',Pct0, ...
%     'ProcessNoise',Qct, ...
%     'MeasurementNoise',R6, ...
%     'StateTransitionFcn',@constturn, ...
%     'MeasurementFcn',measFcnCT);
% 
% %% ===== IMM
% models     = {kfCV, kfCA, kfCT};
% modelNames = {'constvel','constacc','constturn'};
% 
% TPM = [0.997  0.002  0.001;   % CV -> CV strongly
%        0.010  0.988  0.002;   % CA -> stay CA unless clear
%        0.015  0.010  0.975];  % CT -> stay CT a bit, but can return
% 
% initProbs = [0.97; 0.02; 0.01];
% 
% imm = trackingBIMM( ...
%     models, @switchimm, TPM, ...
%     'ModelNames',        modelNames, ...
%     'ModelProbabilities',initProbs, ...
%     'UseWinnerTakesAll', true, ...
%     'YawRateThresh',     0.35, ...   % need a stronger turn to flip to CT
%     'AccelThresh',       3.0,  ...   % need stronger accel to flip to CA
%     'Bias',              0.40, ...   % push prior a bit more
%     'Hysteresis',        0.25, ...   % resist near-ties
%     'StickyDecay',       0.75);      % remember last winner longer
% 
% 
% 
% setMeasurementSizes(imm,6,6);   % detections are 6-D
% end





function [pRow, vRow] = posvelFromTrackCompat(tr)
    x = tr.State(:);
    n = numel(x);
    if n == 6
        % CV
        pRow = [x(1) x(3) x(5)];
        vRow = [x(2) x(4) x(6)];
    elseif n == 7
        % CT (3D)
        pRow = [x(1) x(3) x(6)];
        vRow = [x(2) x(4) x(7)];
    else
        % Fallback to selectors if available (older API)
        try
            pRow = getTrackPositions(tr, [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]);
        catch
            pRow = [NaN NaN NaN];
        end
        try
            vRow = getTrackVelocities(tr, [0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1]);
        catch
            vRow = [0 0 0];
        end
        pRow = pRow(:).';
        vRow = vRow(:).';
    end
end

function Hs = sanitizeHistoryForPlot(H, nearZero, maxJump)
% Replaces bad points with NaNs and breaks before teleports.
% H : Nx2 history [x y] (may already contain NaNs)
% nearZero : coordinate threshold (e.g. 0.2 meters)
% maxJump  : per-step max distance (e.g. 20 meters)

if nargin < 2, nearZero = 0.2; end
if nargin < 3, maxJump  = 20;  end

Hs = H;

% 1) Any point with either coord near zero -> NaN (don’t draw it)
bad0 = (abs(Hs(:,1)) <= nearZero) | (abs(Hs(:,2)) <= nearZero);
Hs(bad0,:) = NaN;

% 2) Break before teleports (distance jump)
idx = find(all(isfinite(Hs),2));          % rows with finite points
if numel(idx) >= 2
    prev = idx(1:end-1);
    next = idx(2:end);
    d = hypot(Hs(next,1)-Hs(prev,1), Hs(next,2)-Hs(prev,2));
    big = d > maxJump;
    % insert NaN right before the big jump
    Hs(prev(big),:) = NaN;
end
end
