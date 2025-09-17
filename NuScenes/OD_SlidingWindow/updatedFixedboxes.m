function [old_boxes, old_label] = updatedFixedboxes( ...
    first_satisfied_index, cluster_history_static, persistence_threshold, ...
    static_box_center, static_bounding_boxes, static_boxe_label, threshold)

    old_boxes = [];
    old_label = [];

    % --- persistence for new labels (internal) ---
    persistent label_counts
    if isempty(label_counts), label_counts = containers.Map('KeyType','char','ValueType','double'); end

    % If no prior frame to compare to, only update counters and exit
    if isempty(first_satisfied_index) || first_satisfied_index < 1 || ...
       first_satisfied_index > numel(cluster_history_static) || ...
       ~isfield(cluster_history_static{first_satisfied_index}, 'static_box_center')
        % Still bump counts for current labels so next frame can pass the gate
        currLabs = static_boxe_label(:);
        for u = 1:numel(currLabs)
            k = num2str(currLabs(u));
            if isKey(label_counts,k)
                label_counts(k) = label_counts(k) + 1;
            else
                label_counts(k) = 1;
            end
        end
        return;
    end

    % Prior frame “accepted” snapshot
    old_boxes = cluster_history_static{first_satisfied_index}.bounding_boxes;
    old_label = cluster_history_static{first_satisfied_index}.static_boxe_label;
    if isempty(old_label), old_label = zeros(0,1); end

    % Ensure column vector
    static_boxe_label = static_boxe_label(:);

    % Ensure old_label has persistence col for compatibility
    if size(old_label,2) < 2
        old_label = [old_label, repmat(persistence_threshold, size(old_label,1), 1)];
    end

    previous_center = cluster_history_static{first_satisfied_index}.static_box_center;
    if isempty(previous_center)
        previous_center = zeros(0,3);
    end

    num_static   = size(static_box_center, 1);
    unmatched_indices = [];

    % Nearest-neighbour gating to previous centers
    for k = 1:num_static
        curr = static_box_center(k, :);
        if isempty(previous_center)
            unmatched_indices(end+1) = k; %#ok<AGROW>
            continue;
        end
        d = sqrt(sum((previous_center(:,1:3) - curr(1:3)).^2, 2));
        mind = min(d);
        if mind >= threshold
            unmatched_indices(end+1) = k; %#ok<AGROW>
        end
    end

    % --- Update per-label persistence counters ---
    % Increment counts for all labels present this frame
    for u = 1:numel(static_boxe_label)
        k = num2str(static_boxe_label(u));
        if isKey(label_counts,k)
            label_counts(k) = label_counts(k) + 1;
        else
            label_counts(k) = 1;
        end
    end
    % (Optional) decay: if a label was previously seen but not present now, gently decay
    prevKeys = keys(label_counts);
    currSet  = string(static_boxe_label(:));
    for kk = 1:numel(prevKeys)
        if ~any(currSet == prevKeys{kk})
            label_counts(prevKeys{kk}) = max(0, label_counts(prevKeys{kk}) - 1);
        end
    end

    % Handle unmatched: only add if label has persisted enough AND not already in old_label
    if ~isempty(unmatched_indices)
        um_boxes  = static_bounding_boxes(unmatched_indices, :);
        um_labels = static_boxe_label(unmatched_indices);   % Nx1

        % keep only labels not already accepted
        already = ismember(um_labels, old_label(:,1));
        um_boxes  = um_boxes(~already, :);
        um_labels = um_labels(~already, :);

        for ii = 1:numel(um_labels)
            k = num2str(um_labels(ii));
            cnt = 0; if isKey(label_counts,k), cnt = label_counts(k); end
            if cnt >= persistence_threshold
                old_boxes = [old_boxes; um_boxes(ii,:)] ; %#ok<AGROW>
                old_label = [old_label; [um_labels(ii), persistence_threshold]]; %#ok<AGROW>
            end
        end
    end

    % Update persistence for accepted labels; remove if they fall below threshold
    for i = size(old_label,1):-1:1
        present_now = ismember(old_label(i,1), static_boxe_label(:));
        if ~present_now
            old_label(i,2) = old_label(i,2) - 1;
            if old_label(i,2) < persistence_threshold
                old_label(i,:) = [];
                old_boxes(i,:) = [];
            end
        else
            old_label(i,2) = persistence_threshold; % reset when present
        end
    end
end
