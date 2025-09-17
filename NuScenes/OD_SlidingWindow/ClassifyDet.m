% function [staticmeas, dynamicmeas,statictime,dynamictime] = ClassifyDet(detections)
% 
% % Initialize output cell arrays outside the loop
% % dynamicmeas = [];
% % staticmeas = [];
% Static = {};
% Dynamic = {};
% numDetections = length(detections);
% 
% for i = 1:numDetections
%     % Extract current measurement
%     meas = detections{i}.Measurement;
%     vx = abs(meas(4));
%     vy = abs(meas(5));
%     vz = abs(meas(6));
% 
%     % Calculate velocity magnitude
%     velocityMagnitude = norm([vx, vy, vz]);
%     % disp(velocityMagnitude)
%     % Classification logic
%     if velocityMagnitude < 0.3  % Static object
%         Static{end+1} = detections{i};
%     else % Dynamic object
%         Dynamic{end+1} = detections{i};
%     end
% end
% if ~isempty(Static) && iscell(Static)
%     staticallDets = [Static{:}];
%     staticmeas = vertcat(staticallDets.Measurement);
%     statictime = vertcat(staticallDets.Time);
% else 
%     staticmeas = [];
%     statictime = [];
% end
% if ~isempty(Dynamic) && iscell(Dynamic)
%     dynamicallDets = [Dynamic{:}];
%     dynamicmeas = vertcat(dynamicallDets.Measurement);
%     dynamictime = vertcat(dynamicallDets.Time);
% else
%     dynamicmeas = [];
%     dynamictime = [];
% end
% end


function [staticmeas, dynamicmeas, statictime, dynamictime] = ClassifyDet(detections)
% detections: cell array of objectDetection

    Static = {}; Dynamic = {};
    numDetections = numel(detections);

    for i = 1:numDetections
        meas = detections{i}.Measurement;
        vx = abs(meas(4)); vy = abs(meas(5)); vz = abs(meas(6));
        speed = norm([vx, vy, vz]);

        if speed < 0.6
            Static{end+1} = detections{i}; %#ok<AGROW>
        else
            Dynamic{end+1} = detections{i}; %#ok<AGROW>
        end
    end

    if ~isempty(Static)
        sAll = [Static{:}];
        staticmeas = double(vertcat(sAll.Measurement));
        statictime = double(vertcat(sAll.Time));   % seconds
    else
        staticmeas = [];
        statictime = [];
    end

    if ~isempty(Dynamic)
        dAll = [Dynamic{:}];
        dynamicmeas = double(vertcat(dAll.Measurement));
        dynamictime = vertcat(dAll.Time);  % seconds
    else
        dynamicmeas = [];
        dynamictime = [];
    end
end
