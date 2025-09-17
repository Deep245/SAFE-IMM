function old_boxes = mergeBoxes(old_boxes, tx, ty, tz)
    if isempty(old_boxes), return; end

    merged_boxes   = [];
    used           = false(size(old_boxes,1),1);
    centers        = [(old_boxes(:,1)+old_boxes(:,4))/2, ...
                      (old_boxes(:,2)+old_boxes(:,5))/2, ...
                      (old_boxes(:,3)+old_boxes(:,6))/2];

    for i = 1:size(old_boxes,1)
        if used(i), continue; end
        group = i;
        ci = centers(i,:);
        for j = i+1:size(old_boxes,1)
            if used(j), continue; end
            cj = centers(j,:);
            if abs(ci(1)-cj(1)) < tx && abs(ci(2)-cj(2)) < ty && abs(ci(3)-cj(3)) < tz
                group(end+1) = j; %#ok<AGROW>
            end
        end
        used(group) = true;
        b = old_boxes(group,:);
        merged_boxes(end+1,:) = [min(b(:,1)), min(b(:,2)), min(b(:,3)), ...
                                 max(b(:,4)), max(b(:,5)), max(b(:,6))]; %#ok<AGROW>
    end

    old_boxes = merged_boxes;
end
