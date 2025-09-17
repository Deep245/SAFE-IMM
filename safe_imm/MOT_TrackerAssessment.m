%% MOT IMM assessment (3 truths) — per-target XY only (truth + all trackers)
% - Per-target "assessment IMM" for diagnostics (ε, WTA, GLR, probs, etc.)
% - One base trackerGNN (linear-SAFE_IMM) used in main loop; two more (nonlinear, baseline)
%   are run after to overlay their outputs on per-target XY plots.
% - RMSE (x y z, vx vy vz) and OSPA(p=1) from tracker outputs (reported).
% - Plots: ONLY per-target XY plots (Truth + Linear + Nonlinear + Baseline)
%          + single combined XY.

clc; clear; rng(7);

%% --------------------------- Timing / sim grid ---------------------------
dt         = 0.1;                  % s
num_frames = 300;
tvec       = (0:num_frames-1)'*dt;
num_targets= 3;

%% --------------------------- Measurement model --------------------------
initFcn = @initIMM_CV_CA_CT_linear;

%posStd = 0.30; velStd = 8.0;
% Measurement model (pos+vel, m=6)
posStd = 2.0; velStd = 0.01; % profile 1
% posStd = 0.01; velStd = 2.0;     % current choice
%posStd = 0.3; velStd = 0.1;  
%posStd = 0.1; velStd = 1.0;  
R_meas = diag([posStd^2 posStd^2 posStd^2 velStd^2 velStd^2 velStd^2]);
mp = []; %#ok<NASGU>

coastFrames = 5;
trk = trackerGNN('FilterInitializationFcn', initFcn, ...
    'AssignmentThreshold', 30, ...
    'ConfirmationThreshold', [5 6], ...
    'DeletionThreshold', coastFrames);

%% --------------------------- Truth generation (3) ------------------------
% State: [x y z vx vy vz]'
Xtruth = zeros(num_frames, 6, num_targets);

% --- Target 1: smooth car-like random (OU) ---
v_min=1.0; v_max=4.0; a_max=3.0; omega_max=0.6;
tau_a=2.0; tau_omega=2.5; sigma_a=1.5; sigma_om=0.5;
bx=120; by=80; k_center=0.003;

x=-60; y=0; psi=0; v=2.5; a=0; om=0;
x1=zeros(num_frames,1); y1=x1; vx1=x1; vy1=x1;
for k=1:num_frames
    a  = a  + (-a/tau_a)*dt     + sigma_a*sqrt(dt)*randn;
    om = om + (-om/tau_omega)*dt+ sigma_om*sqrt(dt)*randn;
    if abs(x) > bx*0.7 || abs(y) > by*0.7
        ang_to_center = atan2(-y, -x);
        dpsi = wrapToPiLocal(ang_to_center - psi);
        om   = om + k_center * dpsi;
    end
    a  = max(min(a,  a_max), -a_max);
    om = max(min(om, omega_max), -omega_max);
    v  = min(max(v + a*dt, v_min), v_max);
    psi= psi + om*dt;
    x  = x + v*cos(psi)*dt; y = y + v*sin(psi)*dt;
    x1(k)=x; y1(k)=y; vx1(k)=v*cos(psi); vy1(k)=v*sin(psi);
end
Xtruth(:,1,1)=x1; Xtruth(:,2,1)=y1; Xtruth(:,3,1)=0;
Xtruth(:,4,1)=vx1;Xtruth(:,5,1)=vy1;Xtruth(:,6,1)=0;

% --- Target 2: straight + loop (two cubic Béziers) ---
x_start=-60; x_join1=-5; x_join2=30; y_peak=18;
h_out=12; h_in=12; top_span_x=6; t_turn_start=12;
idx1 = find(tvec <= t_turn_start); idx2 = setdiff(1:num_frames, idx1);
N1=numel(idx1); N2=numel(idx2);
xA = linspace(x_start, x_join1, max(N1,2)).'; yA = zeros(size(xA));
P0=[x_join1,0]; P1=[x_join1+h_out,0];
P2=[(x_join1+x_join2)/2 - top_span_x, y_peak];
P3=[(x_join1+x_join2)/2,               y_peak];
P4=[(x_join1+x_join2)/2 + top_span_x, y_peak];
P5=[x_join2 - h_in, 0]; P6=[x_join2,0];
N2a=max(floor(N2/2),2); N2b=N2-(N2a-1);
bez = @(P0,P1,P2,P3,t)(1-t).^3.*P0 + 3*(1-t).^2.*t.*P1 + 3*(1-t).*t.^2.*P2 + t.^3.*P3;
tA = linspace(0,1,N2a).'; PA=[bez(P0(1),P1(1),P2(1),P3(1),tA), bez(P0(2),P1(2),P2(2),P3(2),tA)];
tB = linspace(0,1,N2b).'; PB=[bez(P3(1),P4(1),P5(1),P6(1),tB), bez(P3(2),P4(2),P5(2),P6(2),tB)];
PB = PB(2:end,:);
xB=[PA(:,1); PB(:,1)]; yB=[PA(:,2); PB(:,2)];
if N2<=0, xB=[]; yB=[]; end
x2=[xA; xB]; y2=[yA; yB];
if numel(x2) ~= num_frames
    s  = linspace(0,1,numel(x2)).';
    si = linspace(0,1,num_frames).';
    x2 = interp1(s,x2,si,'pchip'); y2 = interp1(s,y2,si,'pchip');
end
vx2 = gradient(x2, dt); vy2 = gradient(y2, dt);
Xtruth(:,1,2)=x2; Xtruth(:,2,2)=y2; Xtruth(:,3,2)=0;
Xtruth(:,4,2)=vx2;Xtruth(:,5,2)=vy2;Xtruth(:,6,2)=0;

% --- Target 3: spline waypoints ---
t_wp = [tvec(1), 6, 12, 18, 24, tvec(end)];
x_wp = [-40, -10, 10, 20, 10, 40]; y_wp = [-20, 10, 25, 0, -15, 10];
ppx = spline(t_wp, x_wp); ppy = spline(t_wp, y_wp);
x3  = ppval(ppx, tvec);   y3  = ppval(ppy, tvec);
vx3 = gradient(x3, dt);   vy3 = gradient(y3, dt);
Xtruth(:,1,3)=x3; Xtruth(:,2,3)=y3; Xtruth(:,3,3)=0;
Xtruth(:,4,3)=vx3;Xtruth(:,5,3)=vy3;Xtruth(:,6,3)=0;

%% --------------------------- Build detections ----------------------------
measCells = cell(num_frames,1);
Z = zeros(6, num_frames, num_targets);
for k=1:num_frames
    dets = cell(1,num_targets);
    for n=1:num_targets
        truth = Xtruth(k,:,n)';
        z = truth + [posStd*randn(3,1); velStd*randn(3,1)];

        % z_clean = truth + [posStd*randn(3,1); velStd*randn(3,1)];
        % 
        % % ---- minimal outlier injection (affects ~10% frames randomly) ----
        % if rand < 0.20
        %     % 6-sigma spike on one random POSITION axis (x/y/z)
        %     j = randi(3);
        %     spike = 6*posStd*randn;
        %     z_clean(j) = z_clean(j) + spike;
        % end
        % 
        % z = z_clean;


        Z(:,k,n) = z;
        dets{n} = objectDetection(tvec(k), z, 'MeasurementNoise', R_meas);
    end
    measCells{k} = dets;
end

%% --------------------------- Assessment IMMs (per target) ----------------
IMMs = cell(1,num_targets);
for n=1:num_targets
    det0 = objectDetection(tvec(1), Z(:,1,n), 'MeasurementNoise', R_meas);
    imm = initFcn(det0);
    if isstruct(imm)
        if isfield(imm,'Filter'), imm = imm.Filter;
        elseif isfield(imm,'imm'), imm = imm.imm;
        else, error('initFcn returned struct without Filter/imm'); end
    end
    try, imm.MeasurementNoise = R_meas; end
    try, imm.pDBG_Enable = true;       end
    IMMs{n} = imm;
end

% Confirm model layouts once (from target 1)
filters = IMMs{1}.TrackingFilters;
dims = cellfun(@(f) numel(f.State), filters(:));
fprintf('Model state dims: %s\n', mat2str(dims));
if any(dims==6) && any(dims==9)
    fprintf('Mapping OK: CV[6]=[x vx y vy z vz], CA[9]=[x vx ax y vy ay z vz az]\n');
else
    fprintf('Check model layouts: expected 6 & 9 states for CV/CA.\n');
end
numModels = numel(filters);

%% --------------------------- Logs (per target) ---------------------------
WinnerIdx  = cell(1,num_targets);
Beps       = cell(1,num_targets);
EpsLimit   = cell(1,num_targets);
WTAflags   = cell(1,num_targets);
GLRsum     = cell(1,num_targets);
Pmax_log   = cell(1,num_targets);
Tail_log   = cell(1,num_targets);
trPbar_log = cell(1,num_targets);
avgd2_log  = cell(1,num_targets);
ModelProbs = cell(1,num_targets);
XestPosIMM = cell(1,num_targets);
XestVelIMM = cell(1,num_targets);

for n=1:num_targets
    WinnerIdx{n}  = zeros(num_frames,1,'uint8');
    Beps{n}       = zeros(num_frames,1,'single');
    EpsLimit{n}   = zeros(num_frames,1,'single');
    WTAflags{n}   = false(num_frames,1);
    GLRsum{n}     = zeros(num_frames,1,'single');
    Pmax_log{n}   = zeros(num_frames,1,'single');
    Tail_log{n}   = zeros(num_frames,1,'single');
    trPbar_log{n} = zeros(num_frames,1,'single');
    avgd2_log{n}  = zeros(num_frames,1,'single');
    ModelProbs{n} = nan(num_frames, numModels, 'single');
    XestPosIMM{n} = nan(num_frames,3);
    XestVelIMM{n} = nan(num_frames,3);
end

% For per-target XY: matched tracker estimate (linear) per truth
EstPosTrackByTruth = cell(1,num_targets);
for n=1:num_targets
    EstPosTrackByTruth{n} = nan(num_frames,3);
end

%% --------------------------- RMSE & OSPA accumulators --------------------
gateDist   = 50;  % m (also non-assignment cost)
ospa_c     = 50;  % cutoff
ospa_p     = 1;   % p=1
ospaCalc = trackOSPAMetric("Metric","OSPA","Distance","posabserr", ...
                           "CutoffDistance",ospa_c,"Order",ospa_p);

% Overall (across targets) per-component RMSE (matched assignments)
se_pos_sum_all = zeros(1,3); % [x y z]
se_vel_sum_all = zeros(1,3); % [vx vy vz]
count_pos_all  = 0;
count_vel_all  = 0;

ospa_all  = nan(num_frames,1);
ospa_loc  = nan(num_frames,1);
ospa_card = nan(num_frames,1);

% Per-target Linear metrics (computed in main loop)
se_pos_lin  = zeros(num_targets,3);   % [x y z]
se_vel_lin  = zeros(num_targets,3);   % [vx vy vz]
cnt_pos_lin = zeros(num_targets,1);
cnt_vel_lin = zeros(num_targets,1);

% Per-target OSPA time series (Linear)
ospa_lin_t    = nan(num_frames, num_targets);
ospa_lin_loc  = nan(num_frames, num_targets);
ospa_lin_card = nan(num_frames, num_targets);

%% --------------------------- Run (FPS measured) --------------------------
t_total = tic; t_imm=0; t_gnn=0; t_log=0;

for k=1:num_frames
    % ---- 1) Assessment IMMs (per target) ----
    tic;
    for n=1:num_targets
        imm = IMMs{n};
        predict(imm, dt);
        z_k = Z(:,k,n);
        [xk, ~, pk] = correct(imm, z_k);

        % map state -> pos/vel indices
        nx = numel(xk);
        if     nx >= 9, ix=[1 4 7]; iv=[2 5 8];  % CA
        elseif nx == 7, ix=[1 3 6]; iv=[2 4 7];  % CT (if used)
        else,           ix=[1 3 5]; iv=[2 4 6];  % CV
        end
        XestPosIMM{n}(k,:) = double(xk(ix)).';
        XestVelIMM{n}(k,:) = double(xk(iv)).';

        % public posteriors
        ModelProbs{n}(k,1:numel(pk)) = single(pk(:)).';
        [~, win] = max(pk); WinnerIdx{n}(k) = uint8(win);

        % GLR proxy
        p = max(pk, realmin('single')); [pmax,imax]=max(p); p(imax)=-inf; p2=max(p);
        GLRsum{n}(k) = (isfinite(p2) && p2>0) * log(double(pmax/p2));

        % debug taps
        try, Beps{n}(k)       = single(imm.pDBG_Last_Beps);     end
        try, EpsLimit{n}(k)   = single(imm.pDBG_Last_EpsLimit); end
        try, WTAflags{n}(k)   = logical(imm.pDBG_WTA);          end
        try, Pmax_log{n}(k)   = single(imm.pDBG_Last_pmax);     end
        try, Tail_log{n}(k)   = single(imm.pDBG_TailMass);      end
        try, trPbar_log{n}(k) = single(imm.pDBG_trPbar);        end
        try, avgd2_log{n}(k)  = single(imm.pDBG_avgd2);         end
    end
    t_imm = t_imm + toc;

    % ---- 2) Run trackerGNN (linear) on all detections for this frame ----
    tic;
    dets = measCells{k};
    out = trk(dets, tvec(k));                  % multi-detection call
    t_gnn = t_gnn + toc;

    % ---- 3) Metrics: assignment, RMSE sums, OSPA ------------------------
    tic;
    % Collect 3D estimates from tracks
    P3 = []; V3 = [];
    if ~isempty(out), out = out([out.IsConfirmed]); end
    if ~isempty(out)
        for q=1:numel(out)
            [pRow, vRow] = posvelFromTrackCompat(out(q));
            p3 = pRow(1:min(3,numel(pRow))); if numel(p3)<3, p3(3)=NaN; end
            v3 = vRow(1:min(3,numel(vRow))); if numel(v3)<3, v3(3)=NaN; end
            P3 = [P3; p3]; V3 = [V3; v3]; %#ok<AGROW>
        end
    end

    % Truth arrays for this frame
    GTpos3 = squeeze(Xtruth(k,1:3,:)).';  % N x 3
    GTvel3 = squeeze(Xtruth(k,4:6,:)).';

    % Assign using 3D position distance
    assignments = zeros(0,2);
    if ~isempty(P3)
        C = pdist2(GTpos3, P3, 'euclidean');
        [assignments, ~, ~] = assignDetectionsToTracks(C, gateDist);
    end

    % Accumulate RMSE, store matched track samples for per-target XY plot
    if ~isempty(assignments)
        for r=1:size(assignments,1)
            n = assignments(r,1); j = assignments(r,2);
            dpos = P3(j,:) - GTpos3(n,:);
            if all(isfinite(dpos))
                se_pos_sum_all = se_pos_sum_all + dpos.^2;
                count_pos_all  = count_pos_all + 1;
            end
            if ~isempty(V3) && all(isfinite(V3(j,:))) && all(isfinite(GTvel3(n,:)))
                dvel = V3(j,:) - GTvel3(n,:);
                se_vel_sum_all = se_vel_sum_all + dvel.^2;
                count_vel_all  = count_vel_all + 1;
            end
            EstPosTrackByTruth{n}(k,:) = P3(j,:);
        end
    end

    % OSPA (3D) over all targets
    truths = repmat(struct('Position',[],'Velocity',[],'PlatformID',[]), num_targets, 1);
    for n=1:num_targets
        truths(n).Position = GTpos3(n,:); truths(n).Velocity = GTvel3(n,:); truths(n).PlatformID = n;
    end
    try
        if isempty(out), trackInput = objectTrack.empty; else, trackInput = out; end
        [d1,loc1,card1] = ospaCalc(trackInput, truths);
        ospa_all(k)=d1; ospa_loc(k)=loc1; ospa_card(k)=card1;
    catch
        d1 = ospaCalc(trackInput, truths);
        ospa_all(k)=d1;
    end

    % ---- Per-target Linear RMSE accumulation & per-target OSPA (Linear) ----
    for n = 1:num_targets
        % Was this truth matched?
        if ~isempty(assignments), idx = find(assignments(:,1)==n, 1, 'first'); else, idx = []; end
        if ~isempty(idx)
            j = assignments(idx,2);
            dpos = P3(j,:) - GTpos3(n,:);
            if all(isfinite(dpos))
                se_pos_lin(n,:)  = se_pos_lin(n,:) + dpos.^2;
                cnt_pos_lin(n)   = cnt_pos_lin(n) + 1;
            end
            if ~isempty(V3) && all(isfinite(V3(j,:))) && all(isfinite(GTvel3(n,:)))
                dvel = V3(j,:) - GTvel3(n,:);
                se_vel_lin(n,:) = se_vel_lin(n,:) + dvel.^2;
                cnt_vel_lin(n)  = cnt_vel_lin(n) + 1;
            end
            trackInput_n = out(j);   % matched-only for OSPA
        else
            trackInput_n = objectTrack.empty; % no match -> cardinality only
        end

        truths1 = struct('Position',GTpos3(n,:), 'Velocity',GTvel3(n,:), 'PlatformID',n);
        try
            [d1n,loc1n,card1n] = ospaCalc(trackInput_n, truths1);
            ospa_lin_t(k,n)    = d1n;
            ospa_lin_loc(k,n)  = loc1n;
            ospa_lin_card(k,n) = card1n;
        catch
            d1n = ospaCalc(trackInput_n, truths1);
            ospa_lin_t(k,n) = d1n;
        end
    end

    t_log = t_log + toc;
end

elapsed = toc(t_total);
fps_total          = num_frames / elapsed;
fps_imm_only       = num_frames / t_imm;
fps_no_trackerGNN  = num_frames / (elapsed - t_gnn);

fprintf('Total: %.2f FPS (elapsed %.3fs)\n', fps_total, elapsed);
fprintf('IMM-only: %.2f FPS (IMM time %.3fs)\n', fps_imm_only, t_imm);
fprintf('Pipeline w/o trackerGNN: %.2f FPS\n', fps_no_trackerGNN);
fprintf('trackerGNN time: %.3fs | logging time: %.3fs\n', t_gnn, t_log);
fprintf('RTF (sim_time/elapsed): %.2fx\n', (tvec(end)-tvec(1)) / elapsed);

%% --------------------------- Prints (concise, per target) ----------------
for n=1:num_targets
    fprintf('\n=== Target %d ===\n', n);
    fprintf('Frames: %d | Elapsed: %.3f s | FPS: %.2f\n', num_frames, elapsed, fps_total);
    sw = [false; diff(double(WinnerIdx{n}))~=0];
    fprintf('Model switches: %d | WTA fired: %d frames\n', nnz(sw), nnz(WTAflags{n}));

    disp('=== IMM Diagnostic Summary ===');
    fprintf('Winner indices (first 20 frames):\n');
    disp(WinnerIdx{n}(1:min(20, num_frames))');

    fprintf('\nε-safe bound values B_eps (first 20):\n');
    disp(Beps{n}(1:min(20, num_frames))');
    epsLimitUsed = max(EpsLimit{n}(EpsLimit{n}>0));
    if isempty(epsLimitUsed) || ~isfinite(epsLimitUsed), epsLimitUsed = NaN; end
    fprintf('Max B_eps observed: %.3f\n', max(Beps{n}));
    fprintf('ε-limit used: %.3f\n', epsLimitUsed);

    violIdx = find(Beps{n} > epsLimitUsed);
    if ~isempty(violIdx)
        fprintf('B_eps > ε-limit at frames (first 20):\n');
        disp(violIdx(1:min(20,numel(violIdx)))');
    else
        fprintf('No WTA ε-limit violations detected.\n');
    end

    numWTA = nnz(WTAflags{n});
    fprintf('\nTotal WTA firings: %d\n', numWTA);
    if numWTA > 0
        fprintf('Frames with WTA fired (first 20):\n');
        disp(find(WTAflags{n}, 20)');
    end

    fprintf('\nFinal model probabilities:\n');
    disp(ModelProbs{n}(num_frames, 1:size(ModelProbs{n},2)));

    fprintf('\nGLR proxy: min=%.3f, max=%.3f, final=%.3f\n', min(GLRsum{n}), max(GLRsum{n}), GLRsum{n}(end));

    disp('=== End of Report ===');
    disp('Inspect plots to verify:');
    disp('- ε-safe gate: B_eps stays under limit when WTA fires (stems on plot 3.1)');
    disp('- WTA events align with sharp posteriors / high GLR (plots 2 & 4)');
    disp('- Winner (plot 2) switches where trajectory curvature/accel changes');

    ok_mask   = WTAflags{n} & (Beps{n} <= EpsLimit{n} | isnan(EpsLimit{n}));
    viol_mask = WTAflags{n} & ~(Beps{n} <= EpsLimit{n} | isnan(EpsLimit{n}));
    fprintf('\nε-gate compliance on WTA frames: %d/%d (%.1f%%)\n', ...
        nnz(ok_mask), nnz(WTAflags{n}), 100*nnz(ok_mask)/max(1,nnz(WTAflags{n})));
    if nnz(viol_mask)
        fprintf('WTA fired while B_eps > limit at frames (first 20):\n');
        disp(find(viol_mask,20)');
    else
        fprintf('No WTA ε-gate violations detected.\n');
    end

    dwell = arrayfun(@(j) nnz(WinnerIdx{n}==j), 1:size(ModelProbs{n},2));
    fprintf('Model dwell frames: %s\n', mat2str(dwell));

    firstWTA = find(WTAflags{n},1,'first'); if isempty(firstWTA), firstWTA = NaN; end
    fprintf('First WTA at frame %g (t=%.2fs)\n', firstWTA, (firstWTA-1)*dt);
    epsMargin = (EpsLimit{n} - Beps{n}); wm = epsMargin(WTAflags{n});
    fprintf('ε-margin on WTA frames: median=%.3f, min=%.3f, max=%.3f\n', median(wm,'omitnan'), min(wm), max(wm));
    pm = max(ModelProbs{n}(WTAflags{n},:),[],2);
    fprintf('p_max on WTA frames: median=%.3f, min=%.3f\n', median(pm,'omitnan'), min(pm));
end

%% --------------------------- RMSE & OSPA summary -------------------------
if count_pos_all>0
    RMSE_pos_all = sqrt(se_pos_sum_all / count_pos_all);
else
    RMSE_pos_all = [NaN NaN NaN];
end
if count_vel_all>0
    RMSE_vel_all = sqrt(se_vel_sum_all / count_vel_all);
else
    RMSE_vel_all = [NaN NaN NaN];
end
fprintf('\n=== RMSE (matched tracks only, all targets) — Linear ===\n');
fprintf('Position RMSE [x y z] (m):      [%0.3f %0.3f %0.3f]\n', RMSE_pos_all);
fprintf('Velocity RMSE [vx vy vz] (m/s): [%0.3f %0.3f %0.3f]\n', RMSE_vel_all);

mOSPA  = mean(ospa_all,'omitnan');
mLoc   = mean(ospa_loc,'omitnan');
mCard  = mean(ospa_card,'omitnan');
fprintf('OSPA p=%d: mean=%0.3f | loc=%0.3f | card=%0.3f (cutoff=%0.1f)\n', ospa_p, mOSPA, mLoc, mCard, ospa_c);

fprintf('\n=== Linear tracker — per-target RMSE & OSPA (from MAIN loop) ===\n');
for n = 1:num_targets
    RMSE_pos_lin_n = sqrt(se_pos_lin(n,:) ./ max(cnt_pos_lin(n),1));
    RMSE_vel_lin_n = sqrt(se_vel_lin(n,:) ./ max(cnt_vel_lin(n),1));
    mOSPA_n  = mean(ospa_lin_t(:,n),   'omitnan');
    mLoc_n   = mean(ospa_lin_loc(:,n), 'omitnan');
    mCard_n  = mean(ospa_lin_card(:,n),'omitnan');
    fprintf('T%d: Pos RMSE [x y z]=[%0.3f %0.3f %0.3f]  ', n, RMSE_pos_lin_n(1), RMSE_pos_lin_n(2), RMSE_pos_lin_n(3));
    fprintf('Vel RMSE [vx vy vz]=[%0.3f %0.3f %0.3f]\n',     RMSE_vel_lin_n(1), RMSE_vel_lin_n(2), RMSE_vel_lin_n(3));
    fprintf('    OSPA p=1: mean=%0.3f | loc=%0.3f | card=%0.3f\n', mOSPA_n, mLoc_n, mCard_n);
end

%% ======================= EXTRA: Nonlinear & Baseline trackers ============
% Re-run detections through two more trackers; collect matched XY estimates.
initLinear   = @initIMM_CV_CA_CT_linear;
initNonlin   = getInitFcnOrFallback('initIMM_CV_CA_CT_nonlinear', initLinear);
initBaseline = getInitFcnOrFallback('initIMM_CV_CA_CT_baseline',  initLinear);

trkB = trackerGNN('FilterInitializationFcn', initNonlin, ...
    'AssignmentThreshold', 30, 'ConfirmationThreshold', [5 6], 'DeletionThreshold', coastFrames);
trkC = trackerGNN('FilterInitializationFcn', initBaseline, ...
    'AssignmentThreshold', 30, 'ConfirmationThreshold', [5 6], 'DeletionThreshold', coastFrames);

% Make sure time starts fresh for each tracker
reset(trkB);
reset(trkC);

% Per-tracker matched positions per truth (1=Nonlinear, 2=Baseline)
trkList = {trkB, trkC};
trackerNames = {'Nonlinear','Baseline'};

% Per-target metrics for EXTRA trackers
se_pos_ex  = zeros(numel(trkList), num_targets, 3);
se_vel_ex  = zeros(numel(trkList), num_targets, 3);
cnt_pos_ex = zeros(numel(trkList), num_targets);
cnt_vel_ex = zeros(numel(trkList), num_targets);
ospa_ex_t    = nan(numel(trkList), num_frames, num_targets);
ospa_ex_loc  = nan(numel(trkList), num_frames, num_targets);
ospa_ex_card = nan(numel(trkList), num_frames, num_targets);

% Per-tracker matched positions per truth
EstPosByTruth_All = cell(1, numel(trkList));
for it=1:numel(trkList)
    EstPosByTruth_All{it} = cell(1,num_targets);
    for n=1:num_targets, EstPosByTruth_All{it}{n} = nan(num_frames,3); end
end

% Overall metrics per extra tracker
RMSEpos_All = nan(numel(trkList),3);
RMSEvel_All = nan(numel(trkList),3);
mOSPA_All   = nan(numel(trkList),1);
mLoc_All    = nan(numel(trkList),1);
mCard_All   = nan(numel(trkList),1);

for it = 1:numel(trkList)
    thisTrk = trkList{it};
    se_pos = zeros(1,3); se_vel = zeros(1,3);
    c_pos = 0; c_vel = 0;
    ospa_all_i = nan(num_frames,1); ospa_loc_i = nan(num_frames,1); ospa_card_i = nan(num_frames,1);

    for k = 1:num_frames
        dets = measCells{k};
        outX = thisTrk(dets, tvec(k));
        if ~isempty(outX), outX = outX([outX.IsConfirmed]); end

        % Gather estimates
        P3=[]; V3=[];
        if ~isempty(outX)
            for q=1:numel(outX)
                [pRow,vRow] = posvelFromTrackCompat(outX(q));
                p3 = pRow(1:min(3,numel(pRow))); if numel(p3)<3, p3(3)=NaN; end
                v3 = vRow(1:min(3,numel(vRow))); if numel(v3)<3, v3(3)=NaN; end
                P3=[P3; p3]; V3=[V3; v3]; %#ok<AGROW>
            end
        end

        GTpos3 = squeeze(Xtruth(k,1:3,:)).';
        GTvel3 = squeeze(Xtruth(k,4:6,:)).';

        assignments=zeros(0,2);
        if ~isempty(P3)
            C = pdist2(GTpos3, P3, 'euclidean');
            [assignments,~,~] = assignDetectionsToTracks(C, gateDist);
        end

        if ~isempty(assignments)
            for r=1:size(assignments,1)
                n = assignments(r,1); j = assignments(r,2);
                dpos = P3(j,:) - GTpos3(n,:);
                if all(isfinite(dpos)), se_pos = se_pos + dpos.^2; c_pos = c_pos+1; end
                if ~isempty(V3) && all(isfinite(V3(j,:))) && all(isfinite(GTvel3(n,:)))
                    dvel = V3(j,:) - GTvel3(n,:); se_vel = se_vel + dvel.^2; c_vel = c_vel+1;
                end
                EstPosByTruth_All{it}{n}(k,:) = P3(j,:);
            end
        end

        % Global OSPA for this tracker
        truths = repmat(struct('Position',[],'Velocity',[],'PlatformID',[]), num_targets, 1);
        for n=1:num_targets
            truths(n).Position = GTpos3(n,:); truths(n).Velocity = GTvel3(n,:); truths(n).PlatformID = n;
        end
        try
            if isempty(outX), trackInput = objectTrack.empty; else, trackInput = outX; end
            [d1,loc1,card1] = ospaCalc(trackInput, truths);
            ospa_all_i(k)=d1; ospa_loc_i(k)=loc1; ospa_card_i(k)=card1;
        catch
            d1 = ospaCalc(trackInput, truths);
            ospa_all_i(k)=d1;
        end

        % ---- Per-target RMSE and per-target OSPA (matched-only) ----
        for n = 1:num_targets
            if ~isempty(assignments), idxn = find(assignments(:,1)==n, 1, 'first'); else, idxn = []; end
            if ~isempty(idxn)
                j = assignments(idxn,2);
                dpos = P3(j,:) - GTpos3(n,:);
                if all(isfinite(dpos))
                    se_pos_ex(it,n,:)  = se_pos_ex(it,n,:) + reshape(dpos.^2,[1 1 3]);
                    cnt_pos_ex(it,n)   = cnt_pos_ex(it,n) + 1;
                end
                if ~isempty(V3) && all(isfinite(V3(j,:))) && all(isfinite(GTvel3(n,:)))
                    dvel = V3(j,:) - GTvel3(n,:);
                    se_vel_ex(it,n,:) = se_vel_ex(it,n,:) + reshape(dvel.^2,[1 1 3]);
                    cnt_vel_ex(it,n)  = cnt_vel_ex(it,n) + 1;
                end
                trackInput_n = outX(j);
            else
                trackInput_n = objectTrack.empty;
            end
            truths1 = struct('Position',GTpos3(n,:), 'Velocity',GTvel3(n,:), 'PlatformID',n);
            try
                [d1n,loc1n,card1n] = ospaCalc(trackInput_n, truths1);
                ospa_ex_t(it,k,n)    = d1n;
                ospa_ex_loc(it,k,n)  = loc1n;
                ospa_ex_card(it,k,n) = card1n;
            catch
                d1n = ospaCalc(trackInput_n, truths1);
                ospa_ex_t(it,k,n) = d1n;
            end
        end
    end

    RMSEpos_All(it,:) = sqrt(se_pos / max(c_pos,1));
    RMSEvel_All(it,:) = sqrt(se_vel / max(c_vel,1));
    mOSPA_All(it)   = mean(ospa_all_i,'omitnan');
    mLoc_All(it)    = mean(ospa_loc_i,'omitnan');
    mCard_All(it)   = mean(ospa_card_i,'omitnan');

    fprintf('\n=== %s tracker — overall RMSE (matched) & OSPA ===\n', trackerNames{it});
    fprintf('Position RMSE [x y z] (m):      [%0.3f %0.3f %0.3f]\n', RMSEpos_All(it,:));
    fprintf('Velocity RMSE [vx vy vz] (m/s): [%0.3f %0.3f %0.3f]\n', RMSEvel_All(it,:));
    fprintf('OSPA p=1: mean=%0.3f | loc=%0.3f | card=%0.3f (cutoff=%0.1f)\n', ...
        mOSPA_All(it), mLoc_All(it), mCard_All(it), ospa_c);
end

for it = 1:numel(trkList)
    fprintf('\n=== %s tracker — per-target RMSE & OSPA ===\n', trackerNames{it});
    for n = 1:num_targets
        RMSE_pos = sqrt(squeeze(se_pos_ex(it,n,:))' ./ max(cnt_pos_ex(it,n),1));
        RMSE_vel = sqrt(squeeze(se_vel_ex(it,n,:))' ./ max(cnt_vel_ex(it,n),1));
        mOSPA  = mean(squeeze(ospa_ex_t(it,:,n)),   'omitnan');
        mLoc   = mean(squeeze(ospa_ex_loc(it,:,n)), 'omitnan');
        mCard  = mean(squeeze(ospa_ex_card(it,:,n)),'omitnan');
        fprintf('T%d: Pos RMSE [x y z]=[%0.3f %0.3f %0.3f]  ', n, RMSE_pos(1), RMSE_pos(2), RMSE_pos(3));
        fprintf('Vel RMSE [vx vy vz]=[%0.3f %0.3f %0.3f]\n',     RMSE_vel(1), RMSE_vel(2), RMSE_vel(3));
        fprintf('    OSPA p=1: mean=%0.3f | loc=%0.3f | card=%0.3f\n', mOSPA, mLoc, mCard);
    end
end

%% --------------------------- Plots (PER-TARGET ONLY) ---------------------
t = tvec(:);

% Colors (consistent across all figures)
colLinear    = [0 0.4470 0.7410];    % blue
colNonlinear = [0.8500 0.3250 0.0980]; % red
colBaseline  = [0.4660 0.6740 0.1880]; % green

% A) Per-target XY: Truth + all tracker outputs (linear/nonlinear/baseline)
for n=1:num_targets
    figure('Name',sprintf('Target %d — Truth vs Trackers (XY)',n), 'NumberTitle','off');
    hold on; grid on; axis equal;
    plot(Xtruth(:,1,n), Xtruth(:,2,n), 'k-', 'LineWidth',2, 'DisplayName','Truth');
    % Linear from main pass:
    plot(EstPosTrackByTruth{n}(:,1), EstPosTrackByTruth{n}(:,2), '-', 'LineWidth',1.5, ...
         'Color',colLinear, 'DisplayName','Linear (matched)');
    % Nonlinear/Baseline from EXTRA runs:
    plot(EstPosByTruth_All{1}{n}(:,1), EstPosByTruth_All{1}{n}(:,2), '-',  'LineWidth',1.5, ...
         'Color',colNonlinear, 'DisplayName','Nonlinear (matched)');
    plot(EstPosByTruth_All{2}{n}(:,1), EstPosByTruth_All{2}{n}(:,2), '-.', 'LineWidth',1.5, ...
         'Color',colBaseline, 'DisplayName','Baseline (matched)');

    legend('Location','best'); xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Target %d: Truth + Trackers (XY)', n));
end

% B) Per-target diagnostics (probs, epsilon gate, GLR)
for n=1:num_targets
    % Model probabilities & switching
    figure('Name',sprintf('Target %d — Model probabilities & switching',n));
    subplot(3,1,1); hold on; grid on;
    for j=1:size(ModelProbs{n},2)
        plot(t, ModelProbs{n}(:,j), '-', 'LineWidth',1);
    end
    ylim([0 1]); ylabel('Prob.'); legend(compose('Model %d',1:size(ModelProbs{n},2)),'Location','eastoutside');
    title('Per-step model probabilities');

    subplot(3,1,2); hold on; grid on;
    stairs(t, double(WinnerIdx{n}), 'LineWidth',1);
    yyaxis right; stem(t(WTAflags{n}), ones(nnz(WTAflags{n}),1), '.'); ylim([0 1.2]);
    ylabel('WTA fire'); yyaxis left; ylabel('Winner idx'); xlabel('Time [s]');
    title('Winner index (stairs) + WTA events (stems)');

    subplot(3,1,3); hold on; grid on;
    plot(t, Pmax_log{n}, '-', 'LineWidth',1);
    plot(t, Tail_log{n}, '-', 'LineWidth',1);
    legend('p_{max}','tail mass (=1-p_{max})','Location','best');
    xlabel('Time [s]'); ylabel('Mass'); title('Sharpness / tail proxy');

    % ε-safe gate diagnostic
    figure('Name',sprintf('Target %d — epsilon-safe gate & separation',n));
    subplot(3,1,1); hold on; grid on;
    plot(t, Beps{n}, '-', 'LineWidth',1);
    if any(EpsLimit{n}>0), plot(t, EpsLimit{n}, '--', 'LineWidth',1); legend('B_{\epsilon}','limit','Location','best');
    else, legend('B_{\epsilon}','Location','best'); end
    stem(t(WTAflags{n}), Beps{n}(WTAflags{n}), '.');
    ylabel('Bound'); title('\epsilon-safe WTA bound vs limit'); xlabel('Time [s]');

    subplot(3,1,2); hold on; grid on;
    plot(t, trPbar_log{n}, '-', 'LineWidth',1);
    ylabel('trace(Pbar)'); title('Mix covariance trace proxy');  % <- fixed label

    subplot(3,1,3); hold on; grid on;
    plot(t, avgd2_log{n}, '-', 'LineWidth',1);
    ylabel('\langle d^2 \rangle'); xlabel('Time [s]');
    title('Avg Mahalanobis^2 to winner');

    % GLR proxy
    figure('Name',sprintf('Target %d — GLR proxy',n)); hold on; grid on;
    plot(t, GLRsum{n}, '-', 'LineWidth',1);
    xlabel('Time [s]'); ylabel('log BF'); title('Maneuver evidence (posterior log-BF)');
end

%% ===== Single plot: all truths (black) + tracker outputs (XY) =====
figure('Name','All Truths + Trackers (XY)','NumberTitle','off'); 
hold on; grid on; axis equal;

% Truths (all black, different line styles)
truth_styles = {'-','--','-.'};
hTruth = gobjects(1,num_targets);
for n = 1:num_targets
    hTruth(n) = plot(Xtruth(:,1,n), Xtruth(:,2,n), truth_styles{min(n,numel(truth_styles))}, ...
        'Color','k','LineWidth',2);
end

% Linear tracker outputs
hL = gobjects(1,num_targets);
for n = 1:num_targets
    EP = EstPosTrackByTruth{n};
    hL(n) = plot(EP(:,1), EP(:,2), '-', 'LineWidth',1.5, 'Color',colLinear);
end

% Nonlinear tracker outputs
hN = gobjects(1,num_targets);
for n = 1:num_targets
    EP = EstPosByTruth_All{1}{n};   % {1}=Nonlinear
    hN(n) = plot(EP(:,1), EP(:,2), '-', 'LineWidth',1.5, 'Color',colNonlinear);
end

% Baseline tracker outputs
hB = gobjects(1,num_targets);
for n = 1:num_targets
    EP = EstPosByTruth_All{2}{n};   % {2}=Baseline
    hB(n) = plot(EP(:,1), EP(:,2), '-.', 'LineWidth',1.5, 'Color',colBaseline);
end

legend([hTruth(1), hL(1), hN(1), hB(1)], ...
       {'Truths','Linear','Nonlinear','Baseline'}, ...
       'Location','bestoutside');
xlabel('x (m)'); ylabel('y (m)');
title('All Truths (black) + Tracker Outputs (XY)');

%% --------------------------- Helpers ------------------------------------
function ang = wrapToPiLocal(ang)
ang = mod(ang + pi, 2*pi) - pi;
end

function [pRow, vRow] = posvelFromTrackCompat(tr)
% Robustly extract [x y z] and [vx vy vz] from objectTrack
pRow = [NaN NaN NaN]; vRow = [NaN NaN NaN];
try
    if isprop(tr,'StateParameters') && ~isempty(tr.StateParameters)
        sp = tr.StateParameters;
        if isstruct(sp)
            if isfield(sp,'PositionSelector') && ~isempty(sp.PositionSelector)
                pRow = (sp.PositionSelector*tr.State).';
            end
            if isfield(sp,'VelocitySelector') && ~isempty(sp.VelocitySelector)
                vRow = (sp.VelocitySelector*tr.State).';
            end
        end
    end
catch
end
x = tr.State(:); n = numel(x);
if any(isnan(pRow)) || any(isnan(vRow))
    switch n
        case 6  % CV: [x vx y vy z vz]
            pRow = [x(1) x(3) x(5)]; vRow = [x(2) x(4) x(6)];
        case 7  % CT: [x vx y vy omega z vz]
            pRow = [x(1) x(3) x(6)]; vRow = [x(2) x(4) x(7)];
        case 9  % CA: [x vx ax y vy ay z vz az]
            pRow = [x(1) x(4) x(7)]; vRow = [x(2) x(5) x(8)];
        otherwise
            try, pRow = getTrackPositions(tr,[1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]); catch, pRow=[NaN NaN NaN]; end
            try, vRow = getTrackVelocities(tr,[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1]); catch, vRow=[0 0 0]; end
            pRow = pRow(:).'; vRow = vRow(:).';
    end
end
end

function f = getInitFcnOrFallback(nameStr, fallbackFcn)
% Returns @nameStr if it exists (on path), else fallbackFcn.
if nargin < 2 || isempty(fallbackFcn)
    error('Fallback initFcn must be provided.');
end
try
    if exist(nameStr,'file') || exist(nameStr,'builtin')
        f = str2func(nameStr);
    else
        f = fallbackFcn;
    end
catch
    f = fallbackFcn;
end
end
