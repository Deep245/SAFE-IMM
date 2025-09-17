% % strict MoMA-M3T-style

clc; clear; close all;
run("varInit.m");   % you provide: measurement_noise, mp, clusterer, thresholds, etc.

%% ------------------- Config (Front-sector, strict metrics) -------------------
scene_dir = 'nuscenes_scene_mat_output';
files = dir(fullfile(scene_dir, 'scene-*.mat'));
if isempty(files), error('No scene-*.mat files found in %s', scene_dir); end

% Tracker variants (comment out what you don’t use)
variants = struct( ...
  'name',     {"linear","nonlinear","baseline"}, ...
  'initFcn',  {@initIMM_CV_CA_CT_linear, @initIMM_CV_CA_CT_nonlinear, @initIMM_CV_CA_CT_baseline}, ...
  'posStd',   {0.6, 1.0, 1.0}, ...
  'velStd',   {17.0, 4.0, 4.0} );

% Matching gate for metrics (meters, XY)
distGateMeters = 4.0;   % try 4–5 for radar

% Base operating point (displayed separately from the sweep)
useConfirmedOnly     = true;
minTrackAgeForEval   = 0;

% FRONT-OF-EGO evaluation wedge (applied to BOTH GT & predictions)
maxRangeEval = 70;      % meters
fovRadarDeg  = 70;     % total FoV => ±50°

% Dynamic-only GT (class-agnostic)
dynOnlyGT            = true;
minSpeedGT           = 0.5;   % m/s considered "moving"
treatUnknownAsDynGT  = false;

% STRICT averaging grid and age sweep
recallGrid = 0:0.01:1;
ageSweep   = 0:6;       % min track age (frames), proxy for confidence

fprintf('Found %d scenes. Evaluating %d variants (front-sector strict metrics, age sweep)...\n', numel(files), numel(variants));

%% ------------------- Result struct -------------------
res = struct('Model',[],'AMOTA',[],'AMOTP_m',[],'MOTA',[],'MOTP_m',[], ...
             'MOTAR',[],'RECALL',[],'GT',[],'MT',[],'ML',[],'TP',[],'FP',[],'FN',[],'IDS',[],'FPS_Hz',[]);

for v = 1:numel(variants)
    name   = variants(v).name;
    initFn = variants(v).initFcn;
    posStd = variants(v).posStd;
    velStd = variants(v).velStd;

    % Accumulators across scenes (base metrics)
    total_GT = 0; total_TP = 0; total_FP = 0; total_FN = 0; total_IDS = 0;
    sum_MOTP_dist = 0; n_frames_total = 0; MT_total = 0; ML_total = 0;

    % Accumulators for AMOTA/MOTAR/AMOTP sweep (age-based)
    S = numel(ageSweep);
    GT_s   = zeros(1,S);
    TP_s   = zeros(1,S);
    FP_s   = zeros(1,S);
    FN_s   = zeros(1,S);
    IDS_s  = zeros(1,S);
    sumDist_s = zeros(1,S);

    t_all_start = tic;

    for k = 1:numel(files)
        scene_file = fullfile(files(k).folder, files(k).name);
        Sdata = load(scene_file);
        samples = Sdata.samples;
        num_frames = numel(samples);

        % Tracker for this variant
        R_meas = diag([posStd^2 posStd^2 posStd^2 velStd^2 velStd^2 velStd^2]);
        coastFrames  = 4;
        trk = trackerGNN( ...
            'FilterInitializationFcn', initFn, ...
            'AssignmentThreshold', 35, ...
            'ConfirmationThreshold', [3 4], ...
            'DeletionThreshold', coastFrames);

        trackerInitialized = false;

        gtByFrame  = cell(1, num_frames);  % struct: id (char), pos (1x2)
        trkByFrame = cell(1, num_frames);  % struct: id (double), pos (1x2), age (double)
        t_by_frame = zeros(1, num_frames); % seconds

        lastGTState = containers.Map('KeyType','char','ValueType','any');  % for speed estimation
        trackerTime = 0;
        prev_tstamp = double(samples{1}.timestamp) * 1e-6;

        for i = 1:num_frames
            s = samples{i};

            % ---------- Time ----------
            t = double(s.timestamp) * 1e-6; t_by_frame(i) = t;
            if i == 1, trackerTime = 0; else, trackerTime = trackerTime + max(t - prev_tstamp, 1e-3); end
            prev_tstamp = t;

            % ---------- Transforms ----------
            R_radar  = quat2rotm(double(s.radar_rotation));   % assumes [w x y z]
            t_radar  = double(s.radar_translation(:));
            R_ego    = quat2rotm(double(s.ego_rotation));
            t_ego    = double(s.ego_translation(:));
            R_global = R_ego * R_radar;                       % RADAR->GLOBAL rotation
            t_global = R_ego*t_radar + t_ego;                 % RADAR->GLOBAL translation

            % ---------- Detections pipeline (RADAR only) ----------
            radarPoints = double(s.radar_points(:, 1:6)); % [x y z vx vy vz] in RADAR frame
            detectionsTbl = processDet(radarPoints, t, 0, 120, -20, 20, -1, 1);

            ObjDet = cell(1, height(detectionsTbl));
            for l = 1:height(detectionsTbl)
                meas = double([detectionsTbl.X(l), detectionsTbl.Y(l), detectionsTbl.Z(l), ...
                               detectionsTbl.Vx(l), detectionsTbl.Vy(l), detectionsTbl.Vz(l)]);
                ObjDet{l} = objectDetection(double(detectionsTbl.t(l)), meas, ...
                    "MeasurementNoise", double(measurement_noise), ...
                    "MeasurementParameters", mp);
            end

            [~, dynamicmeas, ~, dynamictime] = ClassifyDet(ObjDet);

            % Clustering (no viz/history)
            [~, ~, ~, chd] = clusteringRadarPoints(i, window_size, [], dynamictime, ...
                                      [], dynamicmeas, t_by_frame(1:i), ...
                                      clusterer, [], [], [], [], ...
                                      cell(1,num_frames), cell(1,num_frames));

            % Build tracker detections from dynamic clusters (GLOBAL)
            dets = {};
            if ~isempty(chd) && ~isempty(chd{i})
                [~,~,pbb_c] = update_dynamicCandidates( ...
                    i, chd, dynamic_threshold, DynamicBoxThreshold, ...
                    merge_threshold_x, merge_threshold_y, merge_threshold_z, frames_dBox);

                pbb_global = to_global_pbb(pbb_c, R_global, t_global);
                v_global_pbb = (R_global * pbb_c(:,4:6).').';
                pbb_global = double([pbb_global(:,1:3), v_global_pbb]);

                dets = cell(1, size(pbb_global,1));
                for j = 1:size(pbb_global,1)
                    meas = pbb_global(j,:).';
                    dets{j} = objectDetection(double(trackerTime), meas, 'MeasurementNoise', double(R_meas));
                end
            end

            % ---------- Tracker step (first call must have ≥1 detection) ----------
            if ~trackerInitialized
                if isempty(dets)
                    tracks_now = [];
                else
                    tracks_now = trk(dets, trackerTime);
                    trackerInitialized = true;
                end
            else
                tracks_now = trk(dets, trackerTime);
            end

            % Select tracks for base eval (front wedge applies AFTER extraction)
            evalTracks = tracks_now;
            if useConfirmedOnly && ~isempty(evalTracks)
                evalTracks = evalTracks([evalTracks.IsConfirmed]);
            end
            if ~isempty(evalTracks)
                ages = [evalTracks.Age];
                evalTracks = evalTracks(ages > minTrackAgeForEval);
            end

            % Tracks → simple struct (id, pos XY, age) in GLOBAL
            Tstruct_all = struct('id',{},'pos',{},'age',{});
            for it = 1:numel(evalTracks)
                posXY = getXYfromTrack(evalTracks(it));
                if any(~isfinite(posXY)), continue; end
                Tstruct_all(end+1).id  = double(evalTracks(it).TrackID); %#ok<AGROW>
                Tstruct_all(end).pos   = double(posXY(:)).';
                Tstruct_all(end).age   = double(evalTracks(it).Age);
            end

            % FRONT WEDGE filter on predictions (front tracks only, in GLOBAL using ego pose)
            Tstruct = struct('id',{},'pos',{},'age',{});
            for it = 1:numel(Tstruct_all)
                if in_front_wedge_ego(Tstruct_all(it).pos(:), R_ego, t_ego, 0, fovRadarDeg, maxRangeEval)
                    Tstruct(end+1) = Tstruct_all(it); %#ok<AGROW>
                end
            end
            trkByFrame{i} = Tstruct;

            % FRONT WEDGE + DYNAMIC-ONLY GT + SURFACE POINT
            Gstruct = struct('id',{},'pos',{});
            for j = 1:length(s.gt_boxes)
                box = s.gt_boxes{j};
                c   = double(box.center(:));   % GLOBAL [x;y;z]
                if ~in_front_wedge_ego(c(1:2), R_ego, t_ego, 0, fovRadarDeg, maxRangeEval)
                    continue;
                end

                isDyn = true;
                if dynOnlyGT
                    spd = NaN;
                    if isfield(box,'velocity') && ~isempty(box.velocity)
                        vv = double(box.velocity(:));
                        if numel(vv)>=2 && all(isfinite(vv(1:2))), spd = norm(vv(1:2)); end
                    end
                    gid = char(box.instance_token);
                    if ~isfinite(spd)
                        if isKey(lastGTState, gid)
                            prev = lastGTState(gid);
                            dt   = t - prev.t;
                            if dt > 0, spd = norm(c(1:2) - prev.pos) / dt; end
                        end
                    end
                    if ~isfinite(spd), isDyn = treatUnknownAsDynGT; else, isDyn = (spd >= minSpeedGT); end
                    lastGTState(gid) = struct('pos', c(1:2), 't', t);
                end
                if ~isDyn, continue; end

                % Radar-friendly "surface point" instead of box center
                p_surf = gt_surface_point_from_ego(box, t_ego);
                Gstruct(end+1).id  = char(box.instance_token); %#ok<AGROW>
                Gstruct(end).pos   = [p_surf(1) p_surf(2)];
            end
            gtByFrame{i} = Gstruct;
        end % frames

        % ---- Base scene metrics (front wedge already applied above) ----
        metrics_scene = evalMetricsHungarian(gtByFrame, trkByFrame, t_by_frame, distGateMeters, -Inf);

        % Accumulate base
        total_GT        = total_GT  + metrics_scene.GT;
        total_TP        = total_TP  + metrics_scene.TP;
        total_FP        = total_FP  + metrics_scene.FP;
        total_FN        = total_FN  + metrics_scene.FN;
        total_IDS       = total_IDS + metrics_scene.IDS;
        sum_MOTP_dist   = sum_MOTP_dist + metrics_scene.sumDist;
        n_frames_total  = n_frames_total + num_frames;

        MT_total        = MT_total + metrics_scene.MT;
        ML_total        = ML_total + metrics_scene.ML;

        % ---- AMOTA/MOTAR/AMOTP sweep per scene over min track age ----
        for si = 1:S
            thrAge = ageSweep(si);
            ms  = evalMetricsHungarian(gtByFrame, trkByFrame, t_by_frame, distGateMeters, thrAge);
            GT_s(si)      = GT_s(si)      + ms.GT;
            TP_s(si)      = TP_s(si)      + ms.TP;
            FP_s(si)      = FP_s(si)      + ms.FP;
            FN_s(si)      = FN_s(si)      + ms.FN;
            IDS_s(si)     = IDS_s(si)     + ms.IDS;
            sumDist_s(si) = sumDist_s(si) + ms.sumDist;
        end
    end % scenes

    elapsed = toc(t_all_start);
    fps_val = (n_frames_total / max(elapsed, eps));

    %% -------- Base metrics (strict protocol) ----------
    TP = total_TP; FP = total_FP; FN = total_FN; GT = total_GT; IDS = total_IDS;
    RECALL  = TP / max(GT,1);
    MOTA    = 1 - (FN + FP + IDS) / max(GT,1);
    MOTP_m  = (sum_MOTP_dist / max(TP,1));
    r = RECALL;
    MOTAR   = max(0, 1 - (IDS + FP + FN - (1 - r)*GT) / max(r*GT,1));

    %% -------- AMOTA / AMOTP via age sweep ----------
    RECALL_s = TP_s ./ max(GT_s,1);
    MOTAR_s  = zeros(size(RECALL_s));
    for si=1:numel(RECALL_s)
        rsi = RECALL_s(si);
        MOTAR_s(si) = max(0, 1 - (IDS_s(si) + FP_s(si) + FN_s(si) - (1 - rsi).*GT_s(si)) ./ max(rsi.*GT_s(si),1));
    end
    MOTP_s = sumDist_s ./ max(TP_s,1);

    % Build upper envelope in recall space (handles duplicate recalls)
    mask = isfinite(RECALL_s) & isfinite(MOTAR_s);
    if ~any(mask)
        AMOTA = NaN; AMOTP_m = NaN;
    else
        r_vec = RECALL_s(mask);
        motar_vec = MOTAR_s(mask);
        motp_vec  = MOTP_s(mask);

        [r_sorted, ord] = sort(r_vec, 'ascend');
        motar_sorted = motar_vec(ord);
        motp_sorted  = motp_vec(ord);

        [r_u, ~, ic] = unique(r_sorted, 'stable');
        motar_u = accumarray(ic, motar_sorted, [], @max);

        % Co-locate MOTP with the envelope MOTAR choice
        motp_u  = zeros(size(r_u));
        for uu = 1:numel(r_u)
            idxs = find(r_sorted==r_u(uu));
            [~,imax] = max(motar_sorted(idxs));
            motp_u(uu) = motp_sorted(idxs(imax));
        end

        % ---- Robust fixed-grid averaging (no interp1 crash) ----
        if isempty(r_u)
            AMOTA = NaN; AMOTP_m = NaN;

        elseif numel(r_u) == 1
            motar_grid = zeros(size(recallGrid));
            motp_grid  = zeros(size(recallGrid));
            motar_grid(recallGrid <= r_u) = motar_u;
            motp_grid(recallGrid  <= r_u) = motp_u;

            % STRICT: zero beyond max recall already enforced by initialization to 0
            AMOTA   = mean(motar_grid, 'omitnan');
            AMOTP_m = mean(motp_grid,  'omitnan');

        else
            % step-wise hold between known recall points
            motar_grid = interp1(r_u, motar_u, recallGrid, 'previous', 'extrap');
            motp_grid  = interp1(r_u, motp_u,  recallGrid, 'previous', 'extrap');

            % STRICT: zero contribution beyond max achieved recall
            motar_grid(recallGrid > max(r_u)) = 0;
            motp_grid(recallGrid  > max(r_u)) = 0;

            AMOTA   = mean(motar_grid, 'omitnan');
            AMOTP_m = mean(motp_grid,  'omitnan');
        end
    end

    % Save results
    res(v).Model   = name;
    res(v).AMOTA   = AMOTA;     res(v).AMOTP_m = AMOTP_m;
    res(v).MOTA    = MOTA;      res(v).MOTP_m  = MOTP_m;
    res(v).MOTAR   = MOTAR;     res(v).RECALL  = RECALL;
    res(v).GT      = GT;        res(v).MT      = MT_total; res(v).ML = ML_total;
    res(v).TP      = TP;        res(v).FP      = FP;       res(v).FN = FN; res(v).IDS = IDS;
    res(v).FPS_Hz  = fps_val;

    % fprintf(['Variant %-10s | AMOTA %.3f | AMOTP %.2fm | MOTA %.3f | MOTP %.2fm | MOTAR %.3f | ' ...
    %          'Recall %.3f | GT %d | MT %d | ML %d | TP %d | FP %d | FN %d | IDS %d | FPS %.1f\n'], ...
    %         name, AMOTA, AMOTP_m, MOTA, MOTP_m, MOTAR, RECALL, GT, MT_total, ML_total, TP, FP, FN, IDS, fps_val);
end

% ------------------- Table -------------------
% ------------------- Table -------------------
T = struct2table(res);
T = movevars(T, {'Model','AMOTA','AMOTP_m','MOTA','MOTP_m','MOTAR','RECALL','GT','MT','ML','TP','FP','FN','IDS','FPS_Hz'});

disp(' ');
disp('=== CAMERA-BASED METHODS (from MoMA-M3T Table 1) ===');
disp('    Method          AMOTA(%)   AMOTP(m)   MOTA(%)   MOTP(m)   MOTAR(%)    MT     ML');
disp('    ------------------------------------------------------------------------------');
disp('    CenterTrack     4.6        1.543      4.3       0.753     23.1        573    5235');
disp('    TraDeS          5.9        1.49       -         -         -           -      -   ');
disp('    PermaTrack      6.6        1.491      6.0       0.724     32.1        652    5065');
disp('    DEFT            17.7       1.564      15.6      0.770     48.4        1951   3232');
disp('    QD-3DT          21.7       1.550      19.8      0.773     56.3        1893   2970');
disp('    Time3D          21.4       1.36       17.3      0.75      -           -      -   ');
disp('    MoMA-M3T        24.2       1.479      21.3      0.713     58.1        1968   3026');
disp('    MoMA-M3T‡       28.5       1.416      24.6      0.695     62.3        2236   2642');

disp(' ');
disp('=== RADAR FRONT-SECTOR (70° FoV, ≤70 m) | DYNAMIC-ONLY | CLASS-AGNOSTIC ===');
disp('=== STRICT AMOTA/AMOTP: fixed recall grid (0:0.01:1) with zero beyond max recall; MOTA/MOTAR include IDS ===');
disp(T);


[~, idxBest] = max([res.AMOTA] + 1e-6*[res.MOTA] - 1e-6*[res.AMOTP_m]);
fprintf('\nBest model by AMOTA (tie-break MOTA, AMOTP): %s\n', res(idxBest).Model);

%% =================== Helpers ===================

function tf = in_front_wedge_ego(xy_global, R_ego, t_ego, yaw_center, fov_deg, max_range)
% Global XY -> ego, then wedge test centered on ego-forward (yaw_center radians).
p_ego = R_ego.'*([xy_global(:); 0] - t_ego(:));  % global -> ego
x = p_ego(1); y = p_ego(2);
if x <= 0, tf = false; return; end               % only in front of ego x-axis
ang = atan2(y, x) - yaw_center;
ang = atan2(sin(ang), cos(ang));                 % wrap to [-pi,pi]
tf  = (hypot(x,y) <= max_range) && (abs(ang) <= deg2rad(fov_deg/2));
end

function posXY = getXYfromTrack(tr)
% Robust XY extraction across common kinematic layouts.
posXY = [NaN NaN];
try
    % Typical [x vx y vy z vz] selector
    posXY = getTrackPositions(tr,[1 0 0 0 0 0; 0 0 1 0 0 0]);
    return;
catch, end
try
    % [x y z vx vy vz]
    posXY = getTrackPositions(tr,[1 0 0 0 0 0; 0 1 0 0 0 0]);
    return;
catch, end
x = tr.State(:);
switch numel(tr.State)
    case 6, posXY = [x(1) x(3)];
    case 7, posXY = [x(1) x(3)];
end
end

function p_surf = gt_surface_point_from_ego(box, t_ego)
% Nearest BEV point on the oriented 2D GT box surface to the ego position.
% Assumes box.center = [x y z], box.size = [w l h], box.rotation = quaternion (w,x,y,z).
c = double(box.center(:));         % [x;y;z] world
whl = double(box.size(:));         % [w;l;h]  (nuScenes: width, length, height)
w = whl(1); l = whl(2);

Rz = quat2rotm(double(box.rotation));  % 3x3
R2 = Rz(1:2,1:2);                      % yaw rotation in BEV

% ego -> box local coords
d_world = double(t_ego(1:2)) - c(1:2);
d_local = R2' * d_world;

% clamp to rectangle (length along x_local, width along y_local)
p_local = [min(max(d_local(1), -l/2), l/2); ...
           min(max(d_local(2), -w/2), w/2)];

% back to world
p_surf = c(1:2) + R2 * p_local;
end

function M = evalMetricsHungarian(gtByFrame, trkByFrame, t_by_frame, gate, minAge)
% CLEAR MOT with per-frame optimal assignment using assignDetectionsToTracks.
% Filter: keep tracks with Age >= minAge (age sweep is our proxy for confidence).
if nargin < 5 || isempty(minAge),   minAge   = -Inf; end
nF = numel(gtByFrame);

gtState = containers.Map('KeyType','char','ValueType','any'); % id -> struct

TP=0; FP=0; FN=0; IDS=0; FRAG=0; sumDist=0;

for f = 1:nF
    G = gtByFrame{f};
    T = trkByFrame{f};

    % apply minAge filter (for sweep)
    if ~isempty(T) && isfield(T, 'age')
        T = T([T.age] >= minAge);
    end

    Ng = numel(G); Nt = numel(T);

    for g = 1:Ng
        gid = char(G(g).id);
        if ~isKey(gtState, gid)
            gtState(gid) = struct('presentFrames',0,'matchedFrames',0,'everMatched',false,'prevMatched',false,'prevTrackID',NaN);
        end
        S = gtState(gid);
        S.presentFrames = S.presentFrames + 1;
        gtState(gid) = S;
    end

    if Ng==0 && Nt==0, continue; end

    % Distances (XY)
    D = inf(Ng, Nt);
    for g = 1:Ng
        pg = G(g).pos;
        for t = 1:Nt
            pt = T(t).pos;
            D(g,t) = hypot(pg(1)-pt(1), pg(2)-pt(2));
        end
    end

    % Cost & gate
    if Ng>0 && Nt>0
        cost = D;
        cost(D > gate) = gate + 1e6;        % effectively disallow beyond gate
        costOfNonAssignment = gate;
        if exist('assignDetectionsToTracks','file') ~= 2
            error('assignDetectionsToTracks not found. Requires Computer Vision Toolbox.');
        end
        [assignment,~,~] = assignDetectionsToTracks(cost, costOfNonAssignment);
    else
        assignment = zeros(0,2);
    end

    matched_g = false(Ng,1); matched_t = false(Nt,1);
    matches = zeros(0,3);
    for k = 1:size(assignment,1)
        g = assignment(k,1); t = assignment(k,2);
        if D(g,t) <= gate
            matches(end+1,:) = [g, t, D(g,t)]; %#ok<AGROW>
            matched_g(g) = true; matched_t(t) = true;
        end
    end
    if ~isempty(matches), sumDist = sumDist + sum(matches(:,3)); end

    TP = TP + size(matches,1);
    FP = FP + sum(~matched_t);
    FN = FN + sum(~matched_g);

    % IDS bookkeeping
    for g = 1:Ng
        gid = char(G(g).id);
        S = gtState(gid);
        idx = [];
        if ~isempty(matches)
            idx = find(matches(:,1)==g, 1);
        end
        if ~isempty(idx)
            tid   = T(matches(idx,2)).id;
            if S.prevMatched && ~isnan(S.prevTrackID) && S.prevTrackID ~= tid
                IDS = IDS + 1;
            end
            S.matchedFrames = S.matchedFrames + 1;
            S.prevMatched = true; S.prevTrackID = tid; S.everMatched = true;
        else
            S.prevMatched = false;
        end
        gtState(gid) = S;
    end
end

% Finalize per-GT summaries
GT_total = 0; MT=0; ML=0;
keys_ = gtState.keys;
for i = 1:numel(keys_)
    S = gtState(keys_{i});
    GT_total = GT_total + S.presentFrames;
    if S.presentFrames > 0
        r = S.matchedFrames / S.presentFrames;
        if r >= 0.8, MT = MT + 1; end
        if r <= 0.2, ML = ML + 1; end
    end
end

M = struct('GT',GT_total,'TP',TP,'FP',FP,'FN',FN,'IDS',IDS,'FRAG',FRAG,'sumDist',sumDist,'MT',MT,'ML',ML);
end
