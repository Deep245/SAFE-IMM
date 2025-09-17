function [cluster_history_dynamic,pbb,pbb_c] = update_dynamicCandidates( ...
    i,cluster_history_dynamic,dynamic_threshold,DynamicBoxThreshold, ...
    merge_threshold_x,merge_threshold_y,merge_threshold_z,frames_dBox)

    % ===== Normalize thresholds (allow scalar or [dx dy dz]) ====
    if isscalar(DynamicBoxThreshold), DynamicBoxThreshold = [1 1 1]*DynamicBoxThreshold; end
    geom_eps = max(1e-3, dynamic_threshold);  % floor for degenerate sizes

    % ===== Pull previous frame artifacts (for stability/persistence) ====
    if i > 1 && isfield(cluster_history_dynamic{i-1}, 'bounding_boxes')
        previous_bounding_Dboxes   = cluster_history_dynamic{i-1}.bounding_boxes;
        previous_bounding_Dlabels  = cluster_history_dynamic{i-1}.dynamic_box_label;
        previous_bounding_Dcenters = cluster_history_dynamic{i-1}.dynamic_box_center; % Nx6
    else
        previous_bounding_Dboxes = []; previous_bounding_Dlabels = []; previous_bounding_Dcenters = [];
    end

    % sliding label window for persistence
    if i > 1 && isfield(cluster_history_dynamic{i-1}, 'recent_labels')
        recent_labels = cluster_history_dynamic{i-1}.recent_labels;
    else
        recent_labels = cell(1, frames_dBox);
    end

    % per-label prototype (locked size)
    if i > 1 && isfield(cluster_history_dynamic{i-1}, 'proto_size_map')
        proto_size_map = cluster_history_dynamic{i-1}.proto_size_map;
    else
        proto_size_map = containers.Map('KeyType','double','ValueType','any');
    end

    % persistent set from previous frame (for hold-1frame behavior)
    if i > 1 && isfield(cluster_history_dynamic{i-1}, 'pbb_persist')
        prev_persist_boxes = cluster_history_dynamic{i-1}.pbb_persist;
        prev_persist_cent  = cluster_history_dynamic{i-1}.pbb_persist_c;
        prev_persist_lab   = cluster_history_dynamic{i-1}.pbb_persist_labels(:);
    else
        prev_persist_boxes = zeros(0,6);
        prev_persist_cent  = zeros(0,6);
        prev_persist_lab   = zeros(0,1);
    end

    % ===== Early exit if no dynamic labels =====
    if ~isfield(cluster_history_dynamic{i}, 'dynamic_labels')
        % show nothing this frame
        pbb   = zeros(0,6);
        pbb_c = zeros(0,6);
        % keep persistence state
        cluster_history_dynamic{i}.recent_labels       = recent_labels;
        cluster_history_dynamic{i}.proto_size_map      = proto_size_map;
        cluster_history_dynamic{i}.pbb_persist         = zeros(0,6);
        cluster_history_dynamic{i}.pbb_persist_c       = zeros(0,6);
        cluster_history_dynamic{i}.pbb_persist_labels  = zeros(0,1);
        return
    end

    % ===== RAW → BOXES =====
    dynamic_points_current = cluster_history_dynamic{i}.dynamic_points;
    dynamic_labels_current = cluster_history_dynamic{i}.dynamic_labels;

    [dynamic_bounding_boxes, dynamic_box_center_raw, dynamic_box_label] = ...
        createBoundingBoxes(dynamic_points_current, dynamic_labels_current); %#ok<ASGLU>

    % --- merge before IDs
    Dboxes_merged = mergeBoxes(dynamic_bounding_boxes, ...
                               merge_threshold_x, merge_threshold_y, merge_threshold_z);

    % Fallback: if merge killed everything, use raw boxes to ensure a box for every cluster
    if isempty(Dboxes_merged)
        use_boxes   = dynamic_bounding_boxes;
        use_centers = dynamic_box_center_raw(:,1:3);
        use_labels  = dynamic_box_label;
    else
        use_boxes   = Dboxes_merged;
        use_centers = [(Dboxes_merged(:,1)+Dboxes_merged(:,4))/2, ...
                       (Dboxes_merged(:,2)+Dboxes_merged(:,5))/2, ...
                       (Dboxes_merged(:,3)+Dboxes_merged(:,6))/2];
        use_labels  = []; % will be assigned below
    end

    % ===== Stable, one-to-one label assignment (mutual-nearest + gating) =====
    % Goal: avoid switching. We map merged centers to original raw centers
    % with a gate and mutual-nearest rule. Then fill leftovers greedily.
    oc = dynamic_box_center_raw(:,1:3);
    if isempty(oc) || isempty(use_centers)
        assigned_labels = zeros(size(use_centers,1),1);
    else
        C = pdist2(use_centers, oc, 'euclidean');
        % rectangular gate by axis deltas to respect anisotropy
        for r = 1:size(use_centers,1)
            dx = abs(oc(:,1) - use_centers(r,1));
            dy = abs(oc(:,2) - use_centers(r,2));
            dz = abs(oc(:,3) - use_centers(r,3));
            bad = dx > DynamicBoxThreshold(1) | dy > DynamicBoxThreshold(2) | dz > DynamicBoxThreshold(3);
            C(r,bad) = inf;
        end
        assigned_labels = zeros(size(use_centers,1),1);
        % mutual-nearest
        [~, r2o] = min(C,[],2); [~, o2r] = min(C,[],1);
        used_o = false(size(oc,1),1);
        for r = 1:size(C,1)
            o = r2o(r);
            if ~isinf(C(r,o)) && o2r(o) == r && ~used_o(o)
                assigned_labels(r) = dynamic_box_label(o);
                used_o(o) = true;
            end
        end
        % fill leftovers greedily but keep uniqueness of raw centers
        left_r = find(assigned_labels==0);
        for u = left_r.'
            [vals, ord] = sort(C(u,:),'ascend');
            lbl = 0;
            for t = 1:numel(ord)
                if ~isinf(vals(t)) && ~used_o(ord(t))
                    lbl = dynamic_box_label(ord(t));
                    used_o(ord(t)) = true; break;
                end
            end
            if lbl == 0
                % last resort: closest without gate
                d = vecnorm(oc - use_centers(u,:), 2, 2);
                [~, jmin] = min(d); lbl = dynamic_box_label(jmin);
                if used_o(jmin), lbl = dynamic_box_label(jmin); end
            end
            assigned_labels(u) = lbl;
        end
    end

    dynamic_bounding_boxes = use_boxes;
    dynamic_box_label      = assigned_labels;

    % ===== centers + velocities Nx6 (vs previous same label) =====
    if ~isempty(previous_bounding_Dlabels)
        prevL = previous_bounding_Dlabels(:);
        prevC = previous_bounding_Dcenters(:,1:3);
    else
        prevL = []; prevC = [];
    end

    dynamic_box_center = zeros(size(use_centers,1),6);
    dynamic_box_center(:,1:3) = use_centers;
    for r = 1:size(use_centers,1)
        lab = dynamic_box_label(r);
        idx = find(prevL == lab, 1, 'first');
        if ~isempty(idx), v = use_centers(r,:) - prevC(idx,:); else, v = [0 0 0]; end
        dynamic_box_center(r,4:6) = v;
    end

    % ===== (A) SHOW A BOX FOR EVERY CURRENT CLUSTER =====
    % pbb/pbb_c will always reflect current frame clusters (no filtering)
    pbb   = dynamic_bounding_boxes;
    pbb_c = dynamic_box_center;

    % ===== (B) BUILD PERSISTENT SET (2 thresholds) =====
    % T1: frame-count window (N-of-N with sliding buffer)
    pos = mod(i-1, frames_dBox) + 1;
    if numel(recent_labels) ~= frames_dBox, recent_labels = cell(1, frames_dBox); end
    recent_labels{pos} = dynamic_box_label(:);

    filled = min(i, frames_dBox);
    W = unique(vertcat(recent_labels{1:filled}));
    counts = zeros(numel(W),1);
    for t = 1:filled
        counts = counts + ismember(W, recent_labels{t});
    end
    confirmed_labels = W(counts == filled); % present in ALL last N frames

    % T2: geometric stability (center/size close to prototype)
    % Initialize prototypes when 1st confirmed
    for r = 1:size(dynamic_bounding_boxes,1)
        lab = dynamic_box_label(r);
        if ismember(lab, confirmed_labels) && ~isKey(proto_size_map, lab)
            sz = dynamic_bounding_boxes(r,4:6) - dynamic_bounding_boxes(r,1:3);
            sz(sz<=0) = geom_eps;
            proto_size_map(lab) = sz;
        end
    end

    % Filter confirmed labels by geometric proximity to their prototypes
    persist_labels = [];
    persist_boxes  = [];
    persist_cent   = [];

    for r = 1:size(dynamic_bounding_boxes,1)
        lab = dynamic_box_label(r);
        if ~ismember(lab, confirmed_labels), continue; end

        % if prototype missing (rare early), accept but create one
        if ~isKey(proto_size_map, lab)
            proto_size_map(lab) = max(geom_eps, dynamic_bounding_boxes(r,4:6)-dynamic_bounding_boxes(r,1:3));
        end
        proto_sz = proto_size_map(lab);

        % geometric check: center within DynamicBoxThreshold (axis-wise),
        % and size within dynamic_threshold per axis
        ctr = dynamic_box_center(r,1:3);
        sz  = dynamic_bounding_boxes(r,4:6) - dynamic_bounding_boxes(r,1:3);

        % center gating: against last frame same label if available
        pass_center = true;
        idx_prev = find(prev_persist_lab == lab, 1, 'first');
        if ~isempty(idx_prev)
            last_ctr = prev_persist_cent(idx_prev,1:3);
            cdiff = abs(ctr - last_ctr);
            pass_center = all(cdiff <= DynamicBoxThreshold);
        end

        sdiff = abs(sz - proto_sz);
        pass_size = all(sdiff <= max(geom_eps, dynamic_threshold));

        if pass_center && pass_size
            % lock box size to prototype (consistent box)
            bx = [ctr - proto_sz/2, ctr + proto_sz/2];
            persist_boxes  = [persist_boxes; bx]; %#ok<AGROW>
            persist_cent   = [persist_cent; dynamic_box_center(r,1:6)]; %#ok<AGROW>
            persist_labels = [persist_labels; lab]; %#ok<AGROW>
        end
    end

    % Optional hold-1-frame for persistent set if temporarily missing
    % (keeps IDs from flickering out for a single frame)
    if isempty(persist_labels) && ~isempty(prev_persist_lab)
        persist_boxes  = prev_persist_boxes;
        persist_cent   = prev_persist_cent;
        persist_labels = prev_persist_lab;
    end

    % ===== Save everything =====
    cluster_history_dynamic{i}.bounding_boxes     = dynamic_bounding_boxes;
    cluster_history_dynamic{i}.dynamic_box_label  = dynamic_box_label;
    cluster_history_dynamic{i}.dynamic_box_center = dynamic_box_center;

    cluster_history_dynamic{i}.recent_labels       = recent_labels;
    cluster_history_dynamic{i}.proto_size_map      = proto_size_map;

    % (B) persistent outputs
    cluster_history_dynamic{i}.pbb_persist         = persist_boxes;
    cluster_history_dynamic{i}.pbb_persist_c       = persist_cent;
    cluster_history_dynamic{i}.pbb_persist_labels  = persist_labels;

    % (A) per-cluster outputs (required function return)
    % already in pbb, pbb_c

    % legacy locals maintained (not used)
    if ~exist('Dboxes', 'var'); Dboxes = []; end  %#ok<NASGU>
end




