function boxes_g = transform_boxes_to_global(boxes_r, R_global, t_global)
% Convert axis-aligned radar-frame AABBs to global-frame centers,
% preserving width/height so plot function can rotate by yaw.
% Returns [xmin_g, ymin_g, zmin_g, xmax_g, ymax_g, zmax_g] in GLOBAL.

if isempty(boxes_r)
    boxes_g = boxes_r; return;
end

n = size(boxes_r,1);
boxes_g = zeros(n,6);

for k = 1:n
    xmin = boxes_r(k,1); ymin = boxes_r(k,2); zmin = boxes_r(k,3);
    xmax = boxes_r(k,4); ymax = boxes_r(k,5); zmax = boxes_r(k,6);

    w = xmax - xmin; h = ymax - ymin; d = zmax - zmin;
    center_r = [(xmin+xmax)/2; (ymin+ymax)/2; (zmin+zmax)/2];

    center_g = R_global * center_r + t_global;

    boxes_g(k,:) = [ ...
        center_g(1)-w/2, center_g(2)-h/2, center_g(3)-d/2, ...
        center_g(1)+w/2, center_g(2)+h/2, center_g(3)+d/2 ];
end
end




% function boxes_global = transform_boxes_to_global(boxes, rot_global, trans_global)
%     % boxes: Nx6 matrix, each row is [xmin, ymin, zmin, xmax, ymax, zmax]
%     % rot_global: 3x3 global rotation matrix
%     % trans_global: 3x1 translation vector
% 
%     % Initialize an array to store the transformed boxes
%     num_boxes = size(boxes, 1);
%     boxes_global = zeros(num_boxes, 6);  % Preallocate the output array
% 
%     % Loop through each box and apply transformation
%     for i = 1:num_boxes
%         box = boxes(i, :);  % Extract the i-th box [xmin, ymin, zmin, xmax, ymax, zmax]
% 
%         % Extract the min/max coordinates for the current box
%         xmin = box(1);
%         ymin = box(2);
%         zmin = box(3);
%         xmax = box(4);
%         ymax = box(5);
%         zmax = box(6);
% 
%         % Define the 2D corners of the bounding box in the XY plane
%         corners = [
%             xmin,xmax, xmax,xmin;
%             ymin, ymax, ymax, ymin
%         ];
% 
%         % Apply the global rotation and translation to the corners (ignoring Z-axis)
%         global_corners = rot_global(1:2, 1:2) * corners + trans_global(1:2);  % 2x4
% 
%         % Extract the global x and y coordinates
%         x_coords = global_corners(1, :);
%         y_coords = global_corners(2, :);
% 
%         % Return the transformed bounding box in global coordinates as 1x6
%         boxes_global(i, :) = [ min(x_coords), min(y_coords), zmin, max(x_coords), max(y_coords), zmax];
%     end
% end
