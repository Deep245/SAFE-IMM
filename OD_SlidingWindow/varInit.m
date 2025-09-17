% All Control Varibales

global cluster_colors_static cluster_colors_dynamic

% === Measurement parameters and covariance (σx=0.5 m, σv=0.2 m/s) ===
mp = struct('Frame', "Rectangular", ...
    'OriginPosition', zeros(1,3), 'OriginVelocity', zeros(1,3), ...
    'Orientation', eye(3), 'HasVelocity', true, 'IsParentToChild', true);

measurement_noise = diag([0.25 0.25 0.25 0.04 0.04 0.04]); % variances


% Initialize maps (outside of any functions)
cluster_colors_static = containers.Map('KeyType', 'double', 'ValueType', 'any');
cluster_colors_dynamic = containers.Map('KeyType', 'double', 'ValueType', 'any');

% For clustering
clusterer = clusterDBSCAN('MinNumPoints',1,'Epsilon',4, ...
    'EnableDisambiguation',false); % Epsilon 0.2 is best
window_size = 2; % Window size to consider previous frames ( use 3 or 4)


% Static bounding boxes
% First frame to record boundingboxes
first_satisfied_index = [];
threshold_static = 0.001; % Tune/chnage to filter the close x y z points
% Define persistence threshold as a variable
persistence_threshold = 3; % Modify this value as needed, this to add new label if persist for given frame
% Merge boxes
% Initialize thresholds for each dimension
merge_threshold_x = 0.01; % Threshold for X dimension
merge_threshold_y = 0.01; % Threshold for Y dimension
merge_threshold_z = 0.01; % Threshold for Z dimension

% Dynamic bouding boxes
all_dynamic_label = []; % Store all labels
D = []; % Store unique boxes
L = []; % Store unique labels
pbb =[]; % Store previous bounding box
pbb_l = [];
% the framed_box is the most relavent tune factor
frames_dBox = 2; % how many frames to fix the size of bouding box, put 1-6
dynamic_threshold = 0.4; % put the between 0.1-0.4
DynamicBoxThreshold = 1.2; % put between 0.5-0.9



% detection
det = {};

% others
cluster_history_static = cell(0,1);
cluster_history_dynamic = cell(0,1);
% others
staticDetection = [];
DynamicDetection = [];
debug_labels = true;   % set false to silence

%tracks
velScale = 0.3;  % tune arrow length
deleteage = 5;
