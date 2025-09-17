function ObjDet = createDet(detections,measurement_noise,mp)
% Create Object Detection Type and populate ObjDet structure
ObjDet = {};
for i = 1:height(detections)
    ObjDet{i} = objectDetection(detections.t(i), ...
        [detections.X(i), detections.Y(i), detections.Z(i), ...
        detections.Vx(i), detections.Vy(i), detections.Vz(i)], ...
        "MeasurementNoise", measurement_noise, "MeasurementParameters", mp);
end
end