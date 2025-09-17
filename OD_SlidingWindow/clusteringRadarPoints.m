% function [idxD,idxS,cluster_history_static,cluster_history_dynamic] = ...
%     clusteringRadarPoints(i,window_size,statictime,dynamictime,static,dynamic,unique_microseconds,...
%     clusterer,idxD,idxS,dynamic_points_window,static_points_window,cluster_history_static,cluster_history_dynamic)
% if ~isempty(static)
%     % Assign static and dynamic measurements
%     staticX = static(:, 1);
%     staticY = static(:, 2);
%     staticZ = static(:, 3);
%     staticVX = static(:, 4);
%     staticVY = static(:, 5);
%     staticVZ = static(:, 6);
% 
%     % Find static points for this frame (time == unique_microseconds(i))
%     static_points = [staticX(statictime == unique_microseconds(i)), ...
%                  staticY(statictime == unique_microseconds(i)), ...
%                  staticZ(statictime == unique_microseconds(i)), ...
%                  staticVX(statictime == unique_microseconds(i)), ...
%                  staticVY(statictime == unique_microseconds(i)), ...
%                  staticVZ(statictime == unique_microseconds(i))];
%     % Collect static points in the window range [i-window_size, i+window_size]
%     for j = max(1, i - window_size):i
%         static_points_window = [static_points_window; 
%             staticX(statictime == unique_microseconds(j)), ...
%             staticY(statictime == unique_microseconds(j)), ...
%             staticZ(statictime == unique_microseconds(j)), ...
%             staticVX(statictime == unique_microseconds(j)), ...
%             staticVY(statictime == unique_microseconds(j)), ...
%             staticVZ(statictime == unique_microseconds(j))];
%     end
%     % static_labels = dbscan(static_points_window, epsilon, minpts);
%     [idxS,static_labels] = clusterer(static_points_window);
%     static_labels = static_labels';
%     cluster_history_static{i}.static_labels = static_labels;
%     cluster_history_static{i}.static_points = static_points_window;
% else
%     cluster_history_static{i,1} = {}; 
% end
% 
% if ~isempty(dynamic)
%     dynamicX = dynamic(:, 1);
%     dynamicY = dynamic(:, 2);
%     dynamicZ = dynamic(:, 3);
%     dynamicVX = dynamic(:, 4);
%     dynamicVY = dynamic(:, 5);
%     dynamicVZ = dynamic(:, 6);
%     % Find dynamic points for this frame (time == unique_microseconds(i))
%     dynamic_points = [dynamicX(dynamictime == unique_microseconds(i)), ...
%                       dynamicY(dynamictime == unique_microseconds(i)), ...
%                       dynamicZ(dynamictime == unique_microseconds(i)), ...
%                       dynamicVX(dynamictime == unique_microseconds(i)), ...
%                       dynamicVY(dynamictime == unique_microseconds(i)), ...
%                       dynamicVZ(dynamictime == unique_microseconds(i))];
%     % Collect dynamic points in the window range [i-window_size, i+window_size]
%     for j = max(1, i - window_size):i
%         dynamic_points_window = [dynamic_points_window; 
%             dynamicX(dynamictime == unique_microseconds(j)), ...
%             dynamicY(dynamictime == unique_microseconds(j)), ...
%             dynamicZ(dynamictime == unique_microseconds(j)), ...
%             dynamicVX(dynamictime == unique_microseconds(j)), ...
%             dynamicVY(dynamictime == unique_microseconds(j)), ...
%             dynamicVZ(dynamictime == unique_microseconds(j))];
%     end
% 
%     % dynamic_labels = dbscan(dynamic_points_window, epsilon, minpts);
%     [idxD,dynamic_labels] = clusterer(dynamic_points_window);
%     dynamic_labels = dynamic_labels';
%     % Store the updated labels and points in the history
%     cluster_history_dynamic{i}.dynamic_labels = dynamic_labels;
%     cluster_history_dynamic{i}.dynamic_points = dynamic_points_window;
% else
%     cluster_history_dynamic{i,1} = {}; 
% end
% 
% end


function [idxD,idxS,cluster_history_static,cluster_history_dynamic] = ...
    clusteringRadarPoints(i,window_size,statictime,dynamictime,static,dynamic,unique_times,...
    clusterer,idxD,idxS,dynamic_points_window,static_points_window,cluster_history_static,cluster_history_dynamic)

time_tol = 1e-6; % seconds tolerance to match timestamps

% ---- Static ----
if ~isempty(static)
    % time mask for current frame
    mask_now = abs(statictime - unique_times(i)) <= time_tol;

    static_points = [static(mask_now,1:6)];
    % collect window [i-window_size ... i]
    static_points_window = [];
    for j = max(1, i - window_size):i
        mask_j = abs(statictime - unique_times(j)) <= time_tol;
        static_points_window = [static_points_window; static(mask_j,1:6)]; %#ok<AGROW>
    end

    if ~isempty(static_points_window)
        [idxS, static_labels] = clusterer(static_points_window);
        cluster_history_static{i}.static_labels = static_labels(:)';  % row
        cluster_history_static{i}.static_points = static_points_window;
    else
        cluster_history_static{i,1} = {};
    end
else
    cluster_history_static{i,1} = {};
end

% ---- Dynamic ----
if ~isempty(dynamic)
    mask_now = abs(dynamictime - unique_times(i)) <= time_tol;

    dynamic_points = [dynamic(mask_now,1:6)];
    dynamic_points_window = [];
    for j = max(1, i - window_size):i
        mask_j = abs(dynamictime - unique_times(j)) <= time_tol;
        dynamic_points_window = [dynamic_points_window; dynamic(mask_j,1:6)]; %#ok<AGROW>
    end

    if ~isempty(dynamic_points_window)
        [idxD, dynamic_labels] = clusterer(dynamic_points_window);
        cluster_history_dynamic{i}.dynamic_labels = dynamic_labels(:)';
        cluster_history_dynamic{i}.dynamic_points = dynamic_points_window;
    else
        cluster_history_dynamic{i,1} = {};
    end
else
    cluster_history_dynamic{i,1} = {};
end
end
