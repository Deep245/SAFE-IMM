function [cluster_history_static, old_boxes] = update_staticCandiate( ...
    i, cluster_history_static, first_satisfied_index, persistence_threshold, threshold_static, ...
    merge_threshold_x, merge_threshold_y, merge_threshold_z)

    old_boxes = [];

    % Guards
    if isempty(cluster_history_static{i}) || ...
       ~isfield(cluster_history_static{i},'static_points') || ...
       ~isfield(cluster_history_static{i},'static_labels')
        return;
    end

    static_points_current = cluster_history_static{i}.static_points;
    static_labels_current = cluster_history_static{i}.static_labels;

    % Previous-frame snapshots (kept for your archive / compatibility)
    if i > 1 && isfield(cluster_history_static{i-1}, 'bounding_boxes')
        if isempty(first_satisfied_index)
            first_satisfied_index = i - 1;
        end
        previous_bounding_Sboxes   = cluster_history_static{i-1}.bounding_boxes;     %#ok<NASGU>
        previous_bounding_Slabels  = cluster_history_static{i-1}.static_boxe_label;  %#ok<NASGU>
        previous_bounding_Scenters = cluster_history_static{i-1}.static_box_center;  %#ok<NASGU>
    else
        previous_bounding_Sboxes   = []; %#ok<NASGU>
        previous_bounding_Slabels  = []; %#ok<NASGU>
        previous_bounding_Scenters = []; %#ok<NASGU>
    end

    % Create boxes for CURRENT frame clusters
    [static_bounding_boxes, static_box_center, static_box_label] = ...
        createBoundingBoxes(static_points_current, static_labels_current);

    % Store current frame snapshot
    cluster_history_static{i}.bounding_boxes    = static_bounding_boxes;
    cluster_history_static{i}.static_boxe_label = static_box_label;
    cluster_history_static{i}.static_box_center = static_box_center;

    % Update fixed boxes with persistence logic (minimal, no IDs)
    [old_boxes, ~] = updatedFixedboxes( ...
        first_satisfied_index, cluster_history_static, persistence_threshold, ...
        static_box_center, static_bounding_boxes, static_box_label, threshold_static);

    % Merge nearby boxes (distance thresholds)
    if ~isempty(old_boxes)
        old_boxes = mergeBoxes(old_boxes, merge_threshold_x, merge_threshold_y, merge_threshold_z);
    end
end
