clc;
clear;

% put 1 to see plots, else 0
seePlots = 1;

% File selection using a dialog box
[filename, filepath] = uigetfile('*.csv', 'Select a CSV file with radar data');
if isequal(filename, 0)
    disp('No file selected. Exiting...');
    return;
end

% Construct the full path of the selected file
full_filepath = fullfile(filepath, filename);
% Load the CSV file with preserved variable names
data = readtable(full_filepath, 'VariableNamingRule', 'preserve');

% Extract radar data columns
x = data.X; 
y = data.Y; 
z = data.Z; 
v = data.("Radial Velocity"); % Radial velocity from the file
t = data.MicroSeconds; % Assuming this column exists
% Get unique microseconds
unique_microseconds = unique(t);
num_frames = numel(unique_microseconds);
% Initialize a structure to hold points per microsecond
frame_points = cell(num_frames, 1);

% Calculate the total playback duration from the last timestamp
total_duration = unique_microseconds(end); % Use the last timestamp in 't'
% Calculate the time intervals between consecutive frames
time_intervals = diff([0; unique_microseconds]);  
scaled_time_intervals = time_intervals * (total_duration / sum(time_intervals));  

%%
clc

% Init files for variables
run("varInit.m");
% Figure varible init file
run("figInit.m");
% Start the timer for animation
animation_start = tic;

% Resource Consuption
profile on;
% To record memorey usage by mode
[memBefore, sysBefore] = memory;
% CPU load
[~, cpuBefore] = system('wmic cpu get LoadPercentage');
% How many cores
numCores = feature('numCores');

% Loop through frames
for i = 1:num_frames

    % Adjust the orientation and height
    [x,y,z] = adjustOrientation(elev_tilt,az_tilt,sensorHeight,x,y,z);
    % points per frame
    frame_points{i} = [x(t == unique_microseconds(i)), ...
                   y(t == unique_microseconds(i)), ...
                   z(t == unique_microseconds(i)), ...
                   v(t == unique_microseconds(i))]; % Include radial velocity
    % Estimate Velocity
    points = frame_points{i};
    if isempty(points)
        continue;
    end    
    dets = estimateVel(points);
    % Process Points RoI
    detections = processDet(dets, unique_microseconds(i), xLowerBound, xHigerBound, yLowerBound, yHigerBound, zLowerBound, zHigerBound);
    % Create Object Detection Type and populate ObjDet structure
    ObjDet = createDet(detections,measurement_noise,mp);
    % Filter Static and Dynamic Detections
    [static,dynamic,statictime,dynamictime] = ClassifyDet(ObjDet);

    % Clustering variables
    run("clusteringInit.m")
    [idxD,idxS,cluster_history_static,cluster_history_dynamic] = ...
        clusteringRadarPoints(i,window_size,num_frames,statictime,dynamictime,static,dynamic,unique_microseconds,...
        clusterer,idxD,idxS,dynamic_points_window,static_points_window,cluster_history_static,cluster_history_dynamic);
    % Only required for plots
    [all_points,all_labels] = allPoint_forPlots(i,cluster_history_static,cluster_history_dynamic);
    
    % Static detection
    if ~isempty(cluster_history_static{i})
        [cluster_history_static,staticDetection] = update_staticCandiate(i,cluster_history_static,first_satisfied_index,persistence_threshold,threshold_static,...
        merge_threshold_x,merge_threshold_y,merge_threshold_z);
    end

    % Dynamic detection
    if ~isempty(cluster_history_dynamic{i})
        [cluster_history_dynamic,DynamicDetection] = update_dynamicCandidates(i,cluster_history_dynamic,dynamic_threshold,DynamicBoxThreshold,...
        merge_threshold_x,merge_threshold_y,merge_threshold_z,all_dynamic_label,D,frames_dBox,pbb,L);
    end

    % See plots
    if seePlots == 1
        showplots(all_points,all_labels,staticDetection,DynamicDetection,static_labels);
    end
    
    % Pause time calculation
    timeCalculation(i,all_points,num_frames,scaled_time_intervals,animation_start);
end
% End timing the animation loop
animation_duration = toc(animation_start);

% see total profile 
profile viewer;

% % Memory usage after running the model
% [memAfter, sysAfter] = memory;
% % Compute memory used by the model
% modelMemUsed = (memAfter.MemUsedMATLAB - memBefore.MemUsedMATLAB) / 1e6;
% disp(['Memory Used by Model: ', num2str(modelMemUsed), ' MB']);
% 
% % CPU usage after model runs
% [~, cpuAfter] = system('wmic cpu get LoadPercentage');
% % Compute CPU Load Increase Due to Model
% cpuLoadIncrease = sum(cpuAfter) - sum(cpuBefore);
% % Compute Time Per Frame
% timePerFrame = animation_duration / num_frames;
% % Compute CPU Load Per Frame#
% cpuLoadPerFrame = cpuLoadIncrease / num_frames;
% % disp(cpuLoadPerFrame)
% cpuCheck(cpuLoadIncrease);

% Display final timing results
disp(['Algorithm duration: ', num2str(total_duration), ' seconds']);
disp(['Actual duration: ', num2str(animation_duration), ' seconds']);
disp(['Difference from expected: ', num2str(animation_duration - total_duration), ' seconds']);