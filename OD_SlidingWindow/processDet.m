% function detections = processDet(det, time_in_seconds, xLowerBound, xHigerBound, yLowerBound, yHigerBound, zLowerBound, zHigerBound)
%     % Initialize an empty table to store filtered detections
%     D = table([], [], [], [], [],[],[],'VariableNames', {'t', 'X', 'Y', 'Z','Vx','Vy','Vz'});
% 
% 
%     % Loop through each detection in the input structure or table
%     for i = 1:height(det)
%         % Extract detection values
%         X = det(i,1);
%         Y = det(i,2);
%         Z = det(i,3);
% 
%         if abs(det(i,4)) ~= 0
%             Vx = det(i,4);
%             Vy = det(i,5);
%             Vz = det(i,6);
%         else
%             Vx = 0;
%             Vy = 0;
%             Vz = 0;
%         end
% 
%         % Check if detection is within bounds
%         if (X > xLowerBound && X < xHigerBound) && ...
%            (Y > yLowerBound && Y < yHigerBound) && ...
%            (Z > zLowerBound && Z < zHigerBound)
%             % Append valid detection to the table
%             newRow = table(time_in_seconds, X, Y, Z, Vx, Vy, Vz, 'VariableNames', {'t', 'X', 'Y', 'Z', 'Vx','Vy','Vz'});
%             D = [D; newRow];
%         end
%     end
% 
%     % Assign output detections as a table
%     detections = D;
% end

function detections = processDet(det, t_sec, xL, xH, yL, yH, zL, zH)
% det: Nx6 [x y z vx vy vz] in radar frame

    % Mask bounds
    inb = (det(:,1) > xL & det(:,1) < xH) & ...
          (det(:,2) > yL & det(:,2) < yH) & ...
          (det(:,3) > zL & det(:,3) < zH);
    D = det(inb, :);

    % Build table
    n = size(D,1);
    detections = table( ...
        repmat(t_sec, n, 1), D(:,1), D(:,2), D(:,3), D(:,4), D(:,5), D(:,6), ...
        'VariableNames', {'t','X','Y','Z','Vx','Vy','Vz'} );
end
