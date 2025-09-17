% =====================================================
% Tracker Comparision for Linear / Nonlinear / IMM
% =====================================================
clc; clear;
run("varInit.m");   % must define window_size, clusterer, thresholds, etc.

scene   = load('nuscenes_scene_mat_output/scene-0103.mat');
samples = scene.samples;

% ---- Run three tracker variants ----
out.linear    = run_tracker_variant(samples, "linear");
out.nonlinear = run_tracker_variant(samples, "nonlinear");
out.imm       = run_tracker_variant(samples, "baseline");

% ---- Compute metrics (distance threshold in meters) ----
thr = 2.0;
metrics.linear    = compute_mot_metrics(out.linear.GT, out.linear.PR, thr);
metrics.nonlinear = compute_mot_metrics(out.nonlinear.GT, out.nonlinear.PR, thr);
metrics.imm       = compute_mot_metrics(out.imm.GT, out.imm.PR, thr);

% ---- Collect into table ----
T = struct2table(struct( ...
    'Tracker',{'Linear','Nonlinear','IMM'}, ...
    'MOTA',[metrics.linear.MOTA, metrics.nonlinear.MOTA, metrics.imm.MOTA], ...
    'MOTP',[metrics.linear.MOTP, metrics.nonlinear.MOTP, metrics.imm.MOTP], ...
    'Recall',[metrics.linear.Recall, metrics.nonlinear.Recall, metrics.imm.Recall], ...
    'IDS',[metrics.linear.IDS, metrics.nonlinear.IDS, metrics.imm.IDS], ...
    'FP',[metrics.linear.FP, metrics.nonlinear.FP, metrics.imm.FP], ...
    'FN',[metrics.linear.FN, metrics.nonlinear.FN, metrics.imm.FN], ...
    'MT',[metrics.linear.MT, metrics.nonlinear.MT, metrics.imm.MT], ...
    'ML',[metrics.linear.ML, metrics.nonlinear.ML, metrics.imm.ML] ...
));
disp('=== Key MOT Metrics (higher MOTA/Recall, lower MOTP/IDS/FP/FN) ===');
disp(T);

% ---- Deltas vs Linear ----
deltaNL = subtract_metrics(metrics.nonlinear, metrics.linear);
deltaIM = subtract_metrics(metrics.imm,       metrics.linear);
disp('--- Nonlinear - Linear ---'); disp(struct2table(deltaNL));
disp('--- IMM - Linear ---');      disp(struct2table(deltaIM));

% =====================================================
% === Function Definitions ============================
% =====================================================

function out = run_tracker_variant(samples, testVariant)
% Wraps your "working super" code, collects GT + Predictions

num_frames = length(samples);

% === Assign consistent color per object (track) ===
all_ids = {};
for i = 1:num_frames
    for j = 1:length(samples{i}.gt_boxes)
        all_ids{end+1} = samples{i}.gt_boxes{j}.instance_token; %#ok<AGROW>
    end
end
unique_ids = unique(all_ids);
colors     = hsv(length(unique_ids)); %#ok<NASU> 
id_map     = containers.Map(unique_ids, mat2cell(colors, ones(1,length(unique_ids)), 3)); %#ok<NASU>

% === Track histories and clustering state ===
cluster_history_static  = cell(1, num_frames);
cluster_history_dynamic = cell(1, num_frames);
idxS = []; idxD = [];
static_points_window  = [];
dynamic_points_window = [];
unique_microseconds = [];

trackerTime = 0;                          
prev_tstamp = double(samples{1}.timestamp) * 1e-6;
if ~exist('trackHist','var') || isempty(trackHist)
    trackHist = containers.Map('KeyType','double','ValueType','any'); %#ok<NASU>
end

% ===== SELECT IMM VARIANT =====
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

posStd = 1.5;  velStd = 10.0;
R = diag([posStd^2 posStd^2 posStd^2 velStd^2 velStd^2 velStd^2]);
coastFrames  = 3;

trk = trackerGNN( ...
    'FilterInitializationFcn', initFcn, ...
    'AssignmentThreshold', 30, ...
    'ConfirmationThreshold', [2 5], ...
    'DeletionThreshold', coastFrames);

% === Collect results ===
GT = cell(1,num_frames);
PR = cell(1,num_frames);

for i = 1:num_frames
    s = samples{i};

    % === Time ===
    t = double(s.timestamp) * 1e-6;
    unique_microseconds = [unique_microseconds; t]; %#ok<AGROW>
    if i == 1
        trackerTime = 0;
    else
        Ts = max(t - prev_tstamp, 1e-3);
        trackerTime = trackerTime + Ts;
    end
    prev_tstamp = t;

    % === Transforms ===
    R_radar = quat2rotm(s.radar_rotation);
    t_radar = s.radar_translation(:);
    R_ego   = quat2rotm(s.ego_rotation);
    t_ego   = s.ego_translation(:);
    R_global = R_ego * R_radar;
    t_global = R_ego * t_radar + t_ego;

    % === Bounds & detection table (x y z vx vy vz) ===
    radarPoints = s.radar_points(:, 1:6);
    detections = processDet(radarPoints, t, 0,80,-20,20,-1,1);

    % Build objectDetections for classification (no external vars)
    ObjDet = cell(1, height(detections));
    for l = 1:height(detections)
        ObjDet{l} = objectDetection( ...
            detections.t(l), ...
            [detections.X(l), detections.Y(l), detections.Z(l), ...
             detections.Vx(l), detections.Vy(l), detections.Vz(l)]);
    end

    % === Static / Dynamic classification ===
    [staticmeas, dynamicmeas, statictime, dynamictime] = ClassifyDet(ObjDet);

    % === Clustering (sliding window over time indices) ===
    [idxD, idxS, cluster_history_static, cluster_history_dynamic] = ...
        clusteringRadarPoints(i, evalin('base','window_size'), statictime, dynamictime, ...
                              staticmeas, dynamicmeas, unique_microseconds, ...
                              evalin('base','clusterer'), idxD, idxS, ...
                              dynamic_points_window, static_points_window, ...
                              cluster_history_static, cluster_history_dynamic);

    % === Dynamic detection → tracker detections ===
    dets = {};
    if ~isempty(cluster_history_dynamic{i})
        [cluster_history_dynamic,DynamicDetection,pbb_c] = update_dynamicCandidates( ...
            i, cluster_history_dynamic, ...
            evalin('base','dynamic_threshold'), evalin('base','DynamicBoxThreshold'), ...
            evalin('base','merge_threshold_x'), evalin('base','merge_threshold_y'), ...
            evalin('base','merge_threshold_z'), evalin('base','frames_dBox'));
        DynamicDetection = transform_boxes_to_global(DynamicDetection, R_global, t_global); %#ok<NASU>

        % velocities rotate only
        v_global_pbb = (R_global * pbb_c(:,4:6).').';
        pbb_global   = to_global_pbb(pbb_c, R_global, t_global);
        pbb_global   = [pbb_global(:,1:3), v_global_pbb];

        dets = cell(1, size(pbb_global,1));
        for j = 1:size(pbb_global,1)
            meas = pbb_global(j,:).';
            dets{j} = objectDetection(trackerTime, meas, 'MeasurementNoise', R);
        end
    end

    % === Step tracker ===
    trac = trk(dets, trackerTime);

    % === Keep confirmed+aged tracks for eval ===
    pr = struct('id', {}, 'c', {});
    if ~isempty(trac)
        trac = trac([trac.IsConfirmed]);
        if ~isempty(trac)
            ages = [trac.Age];
            trac = trac(ages > 3);
            for k = 1:numel(trac)
                st = trac(k).State; % [x;y;z;vx;vy;vz]
                pr(end+1).id = double(trac(k).TrackID); %#ok<AGROW>
                pr(end).c    = st(1:2).';
            end
        end
    end
    PR{i} = pr;

    % === Collect GT with front-sector gates (radar OR camera) ===
    gt = struct('id', {}, 'c', {});
    maxRange = 100; halfFovRadar = deg2rad(30); halfFovCam = deg2rad(30);
    yawRadar = atan2(R_global(2,1), R_global(1,1));
    yawEgo   = atan2(R_ego(2,1), R_ego(1,1));
    if isfield(s, 'camera_rotation')
        R_cam     = quat2rotm(s.camera_rotation);
        R_cam_glo = R_ego * R_cam;
        yawCam    = atan2(R_cam_glo(2,3), R_cam_glo(1,3));
    else
        yawCam = yawEgo;
    end
    for jj = 1:length(s.gt_boxes)
        box   = s.gt_boxes{jj};
        center = box.center(:);
        dx = center(1) - t_ego(1); dy = center(2) - t_ego(2);
        if hypot(dx,dy) > maxRange, continue; end
        theta = atan2(dy, dx);
        inRadar  = abs(wrapToPi(theta - yawRadar)) <= halfFovRadar;
        inCamera = abs(wrapToPi(theta - yawCam))   <= halfFovCam;
        if ~(inRadar || inCamera), continue; end
        gt(end+1).id = string2hash(box.instance_token); %#ok<AGROW>
        gt(end).c    = center(1:2).';
    end
    GT{i} = gt;
end

out.GT = GT; out.PR = PR;
end

% --- Compute MOT metrics (uses assignDetectionsToTracks) ---
function M = compute_mot_metrics(GT, PR, thr)
numF = numel(GT); TP=0;FP=0;FN=0;IDS=0;sumDist=0;numTP=0;
gtIds=[]; for t=1:numF, gtIds=[gtIds,[GT{t}.id]]; end; gtIds=unique(gtIds);

gtSeen = containers.Map('KeyType','double','ValueType','double');
gtCov  = containers.Map('KeyType','double','ValueType','double');
gtLast = containers.Map('KeyType','double','ValueType','double');

for t = 1:numF
    g = GT{t}; p = PR{t}; gN = numel(g); pN = numel(p);

    for i = 1:gN
        gid = g(i).id;
        if ~isKey(gtSeen,gid), gtSeen(gid)=0; end
        gtSeen(gid) = gtSeen(gid)+1;
    end
    if gN==0 && pN==0, continue; end

    if gN==0
        FP = FP + pN;    % all predictions are false positives
        continue;
    elseif pN==0
        FN = FN + gN;    % all GT are missed
        continue;
    end

    % Distance matrix
    D = zeros(gN,pN);
    for i=1:gN, for j=1:pN, D(i,j)=norm(g(i).c - p(j).c); end, end

    % Thresholding + Hungarian via assignDetectionsToTracks
    Dth = D; Dth(Dth>thr) = inf;
    % costOfNonAssignment = thr ensures unassignment if all costs > thr
    [assignments, unTrk, unDet] = assignDetectionsToTracks(Dth, thr); %#ok<ASGLU>

    % True positives and distance
    for a = 1:size(assignments,1)
        gi = assignments(a,1);
        pj = assignments(a,2);
        if isfinite(D(gi,pj)) && D(gi,pj) <= thr
            TP = TP + 1;
            sumDist = sumDist + D(gi,pj);
            numTP   = numTP + 1;

            gid = g(gi).id; pid = p(pj).id;
            if ~isKey(gtCov,gid), gtCov(gid)=0; end
            gtCov(gid) = gtCov(gid)+1;

            if isKey(gtLast,gid) && gtLast(gid)~=pid, IDS = IDS + 1; end
            gtLast(gid) = pid;
        end
    end

    % FP: unmatched predictions
    FP = FP + numel(unDet);
    % FN: unmatched GT
    FN = FN + numel(unTrk);
end

% MT/ML
MT=0; ML=0;
for k=1:numel(gtIds)
    gid = gtIds(k);
    present = 0; covered = 0;
    if isKey(gtSeen,gid), present = gtSeen(gid); end
    if isKey(gtCov, gid), covered = gtCov(gid); end
    if present==0, continue; end
    covr = covered/present;
    if covr>=0.8, MT=MT+1; elseif covr<=0.2, ML=ML+1; end
end

GTtotal = 0; for t=1:numF, GTtotal = GTtotal + numel(GT{t}); end
if GTtotal==0
    MOTA=0; Recall=0;
else
    MOTA = 1 - (FN + FP + IDS)/GTtotal;
    Recall = (TP)/GTtotal;
end
MOTP = (numTP>0) * (sumDist / max(numTP,1));

M = struct('MOTA',MOTA,'MOTP',MOTP,'Recall',Recall, ...
           'IDS',IDS,'FP',FP,'FN',FN,'MT',MT,'ML',ML);
end

% --- Utility: subtract metrics ---
function d=subtract_metrics(a,b)
d=struct('MOTA',a.MOTA-b.MOTA,'MOTP',a.MOTP-b.MOTP,'Recall',a.Recall-b.Recall, ...
         'IDS',a.IDS-b.IDS,'FP',a.FP-b.FP,'FN',a.FN-b.FN,'MT',a.MT-b.MT,'ML',a.ML-b.ML);
end

% --- Deterministic numeric ID from string token (double-only math) ---
function h = string2hash(str)
str = char(str);
h = sum( double(str) .* (1:numel(str)) );
end
