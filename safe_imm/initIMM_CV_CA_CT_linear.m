function imm = initIMM_CV_CA_CT_linear(detection)
% INIT for trackerGNN: returns a single tracking filter object.
% detection.Measurement: [x y z vx vy vz], detection.MeasurementNoise: 6x6



keepCT = false;

m  = double(detection.Measurement(:));
R6 = double(detection.MeasurementNoise);

vx0 = (numel(m)>=4)*m(4);
vy0 = (numel(m)>=5)*m(5);
vz0 = (numel(m)>=6)*m(6);

H6_CV = [ 1 0 0 0 0 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 0 1 ];

H6_CA = [ 1 0 0 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0 0;
          0 0 0 0 0 0 1 0 0;
          0 1 0 0 0 0 0 0 0;
          0 0 0 0 1 0 0 0 0;
          0 0 0 0 0 0 0 1 0];

Fcv = @(dt) blkdiag([1 dt;0 1],[1 dt;0 1],[1 dt;0 1]);
Qcv_axis = @(dt,q) q * [dt^4/4, dt^3/2; dt^3/2, dt^2];
Qcv = @(dt,q) blkdiag(Qcv_axis(dt,q), Qcv_axis(dt,q), Qcv_axis(dt,q));
qCV = 0.5;

xcv0 = [m(1);vx0;  m(2);vy0;  m(3);vz0];
Pcv0 = diag([4 25  4 25  4 25]);

kfCV = trackingKF( ...
    "MotionModel","Custom", ...
    "StateTransitionModel", Fcv(1), ...
    "MeasurementModel",     H6_CV, ...
    "State",                xcv0, ...
    "StateCovariance",      Pcv0, ...
    "ProcessNoise",         Qcv(1,qCV), ...
    "MeasurementNoise",     R6);

F3  = @(dt) [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];
Fca = @(dt) blkdiag(F3(dt),F3(dt),F3(dt));
Qca_axis = @(dt,q) q * [dt^5/20, dt^4/8, dt^3/6; dt^4/8, dt^3/3, dt^2/2; dt^3/6, dt^2/2, dt];
Qca = @(dt,q) blkdiag(Qca_axis(dt,q), Qca_axis(dt,q), Qca_axis(dt,q));
qCA = 0.2;

xca0 = [m(1);vx0;0;  m(2);vy0;0;  m(3);vz0;0];
Pca0 = diag([4 25 25,  4 25 25,  9 36 36]);

kfCA = trackingKF( ...
    "MotionModel","Custom", ...
    "StateTransitionModel", Fca(1), ...
    "MeasurementModel",     H6_CA, ...
    "State",                xca0, ...
    "StateCovariance",      Pca0, ...
    "ProcessNoise",         Qca(1,qCA), ...
    "MeasurementNoise",     R6);

if keepCT
    Hct6 = [ 1 0 0 0 0 0 0
             0 0 1 0 0 0 0
             0 0 0 0 0 1 0
             0 1 0 0 0 0 0
             0 0 0 1 0 0 0
             0 0 0 0 0 0 1 ];
    xct0 = [m(1); vx0; m(2); vy0; 0; m(3); vz0];
    Pct0 = diag([4 25  4 25  (pi/2)^2  9 36]);
    Qct  = blkdiag(0.2*eye(2), 0.2*eye(2), 0.02, 0.4*eye(2));
    ekfCT = trackingEKF( ...
        "State",                   xct0, ...
        "StateCovariance",         Pct0, ...
        "StateTransitionFcn",      @(x,dt) constturn(x,dt), ...
        "ProcessNoise",            Qct, ...
        "MeasurementFcn",          @(x) Hct6*x, ...
        "MeasurementJacobianFcn",  @(x) Hct6, ...
        "MeasurementNoise",        R6);
    models = {kfCV, kfCA, ekfCT};
    names  = {'constvel','constacc','constturn'};
    TPM    = [0.986 0.010 0.004; 0.010 0.986 0.004; 0.015 0.015 0.970];
else
    models = {kfCV, kfCA};
    names  = {'constvel','constacc'};
    TPM    = [0.992 0.008; 0.015 0.985];
end

% IMPORTANT: use the exact class name that matches its file and classdef
imm = trackingDynamicLinearIMM_novel( ...
    models, ...
    @switchimm, TPM, ...
    'YawRateThresh', 0.12, ...
    'AccelThresh',   1.8, ...
    'Bias',          0.22, ...
    'Hysteresis',    0.04, ...
    'StickyDecay',   0.25, ...
    'ModelNames',    names);

%setMeasurementSizes(imm, 6, 6);
end

