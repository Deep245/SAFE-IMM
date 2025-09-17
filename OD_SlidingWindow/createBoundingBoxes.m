function [bounding_boxes, cluster_centers, unique_labels] = createBoundingBoxes(points, labels)
% points: Nx6 [x y z vx vy vz] or Nx3; labels: Nx1 or 1xN

    if isempty(points) || isempty(labels)
        bounding_boxes  = zeros(0,6);
        cluster_centers = zeros(0,3);
        unique_labels   = zeros(0,1);
        return;
    end

    labels = labels(:);
    valid_mask   = labels > 0;
    valid_points = points(valid_mask, :);
    valid_labels = labels(valid_mask);

    if isempty(valid_points)
        bounding_boxes  = zeros(0,6);
        cluster_centers = zeros(0,3);
        unique_labels   = zeros(0,1);
        return;
    end

    unique_labels   = unique(valid_labels);
    nC = numel(unique_labels);

    bounding_boxes  = zeros(nC, 6);
    cluster_centers = zeros(nC, 3);

    % --- minimum box side lengths (meters) ---
    min_box = [1.5, 1.5, 1.5];  % tweak if needed

    for k = 1:nC
        lab = unique_labels(k);
        cp  = valid_points(valid_labels == lab, :);
        xyz = cp(:, 1:min(3,size(cp,2)));

        % Center of the cluster
        c = mean(xyz, 1);
        cluster_centers(k,:) = c;

        if size(xyz,1) == 1
            % Single-point cluster: directly use min_box around the point
            w = min_box(1); h = min_box(2); d = min_box(3);
        else
            % Multi-point cluster: raw extents
            xmin0 = min(xyz(:,1)); xmax0 = max(xyz(:,1));
            ymin0 = min(xyz(:,2)); ymax0 = max(xyz(:,2));
            zmin0 = min(xyz(:,3)); zmax0 = max(xyz(:,3));

            w0 = max(xmax0 - xmin0, 0);
            h0 = max(ymax0 - ymin0, 0);
            d0 = max(zmax0 - zmin0, 0);

            % Clamp each side to at least min_box (around the center)
            w = max(w0, min_box(1));
            h = max(h0, min_box(2));
            d = max(d0, min_box(3));
        end

        % Reconstruct bbox around the center with clamped sizes
        xmin = c(1) - w/2; xmax = c(1) + w/2;
        ymin = c(2) - h/2; ymax = c(2) + h/2;
        zmin = c(3) - d/2; zmax = c(3) + d/2;

        bounding_boxes(k,:) = [xmin, ymin, zmin, xmax, ymax, zmax];
    end
end
