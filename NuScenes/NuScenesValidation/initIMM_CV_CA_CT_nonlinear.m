%--- IMM initializer: CV + CA (+ optional CT), 6D measurements ---
function imm = initIMM_CV_CA_CT_nonlinear(detection)
% Expects a detection with 6-D measurement: [x y z vx vy vz]

    % Measurement + class helper
    R6 = detection.MeasurementNoise;
    m  = detection.Measurement(:);              % [x y z vx vy vz]
    clsLike = m(1);                             % for casting (single/double)
    vx0 = (numel(m)>=4) * m(4);
    vy0 = (numel(m)>=5) * m(5);
    vz0 = (numel(m)>=6) * m(6);

    % Cast helper to match numeric class
    castLike = @(A,x) cast(A,'like',x);

    %% ======= CV (α–β): state [x vx y vy z vz], meas 6D =======
    Fcv_mat = @(dt) [ 1  dt   0   0   0   0
                      0   1   0   0   0   0
                      0   0   1  dt   0   0
                      0   0   0   1   0   0
                      0   0   0   0   1  dt
                      0   0   0   0   0   1 ];
    Hcv6_mat = [ 1 0  0 0  0 0   % x
                 0 0  1 0  0 0   % y
                 0 0  0 0  1 0   % z
                 0 1  0 0  0 0   % vx
                 0 0  0 1  0 0   % vy
                 0 0  0 0  0 1]; % vz
    % Q for dt≈1 (white-noise accel, lifted to 3 axes)
    Qcv_dt1 = @(q) blkdiag(q * ([0.5;1]*[0.5 1]), ...
                           q * ([0.5;1]*[0.5 1]), ...
                           q * ([0.5;1]*[0.5 1]));

    xcv0 = [m(1); vx0;  m(2); vy0;  m(3); vz0];
    Pcv0 = diag([ 4 25  4 25  4 25 ]);
    qCV  = 0.5;  % TUNE

    ekfCV = trackingEKF( ...
        'State',                      xcv0, ...
        'StateCovariance',            Pcv0, ...
        'ProcessNoise',               castLike(Qcv_dt1(qCV), clsLike), ...
        'MeasurementNoise',           R6, ...
        'StateTransitionFcn',         @(x,dt) castLike(Fcv_mat(dt),x)*x, ...
        'StateTransitionJacobianFcn', @(x,dt) castLike(Fcv_mat(dt),x), ...
        'MeasurementFcn',             @(x)    castLike(Hcv6_mat,x)*x, ...
        'MeasurementJacobianFcn',     @(x)    castLike(Hcv6_mat,x) );

    %% ======= CA (α–β–γ): state [x vx ax y vy ay z vz az], meas 6D =======
    F3 = @(dt) [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];
    Fca_mat = @(dt) blkdiag(F3(dt),F3(dt),F3(dt));

    Hca6_mat = [ 1 0 0   0 0 0   0 0 0   % x
                 0 0 0   1 0 0   0 0 0   % y
                 0 0 0   0 0 0   1 0 0   % z
                 0 1 0   0 0 0   0 0 0   % vx
                 0 0 0   0 1 0   0 0 0   % vy
                 0 0 0   0 0 0   0 1 0]; % vz

    % Q for dt≈1 (white jerk in 1D, lifted to 3 axes)
    Q1_wj = [ 1/20  1/8  1/6
              1/8   1/3  1/2
              1/6   1/2  1   ];
    Qca_dt1 = @(q) blkdiag(q*Q1_wj, q*Q1_wj, q*Q1_wj);

    xca0 = [m(1); vx0; 0;   m(2); vy0; 0;   m(3); vz0; 0];
    Pca0 = diag([ 4 25 25,  4 25 25,  9 36 36 ]);
    qCA  = 0.2;  % TUNE

    ekfCA = trackingEKF( ...
        'State',                      xca0, ...
        'StateCovariance',            Pca0, ...
        'ProcessNoise',               castLike(Qca_dt1(qCA), clsLike), ...
        'MeasurementNoise',           R6, ...
        'StateTransitionFcn',         @(x,dt) castLike(Fca_mat(dt),x)*x, ...
        'StateTransitionJacobianFcn', @(x,dt) castLike(Fca_mat(dt),x), ...
        'MeasurementFcn',             @(x)    castLike(Hca6_mat,x)*x, ...
        'MeasurementJacobianFcn',     @(x)    castLike(Hca6_mat,x) );

    %% ======= Optional CT model (lightweight) =======
    keepCT = true;
    models    = {ekfCV, ekfCA};
    names     = ["constvel","constacc"];
    TPM       = [0.992 0.008;
                 0.015 0.985];
    initProbs = [0.90 0.10];

    if keepCT
        % State: [x vx y vy w z vz]  (w = yaw rate)
        Hct6_mat = [ 1 0 0 0 0 0 0
                     0 0 1 0 0 0 0
                     0 0 0 0 0 1 0
                     0 1 0 0 0 0 0
                     0 0 0 1 0 0 0
                     0 0 0 0 0 0 1];
        xct0 = [m(1); vx0; m(2); vy0; 0; m(3); vz0];
        Pct0 = diag([ 4 25  4 25  (pi/2)^2  9 36 ]);

        Qct  = blkdiag( ...
            0.2*eye(2,'like',clsLike), ... % x,vx
            0.2*eye(2,'like',clsLike), ... % y,vy
            cast(0.02,'like',clsLike), ... % w
            0.4*eye(2,'like',clsLike));    % z,vz

        ekfCT = trackingEKF( ...
            'State',                   xct0, ...
            'StateCovariance',         Pct0, ...
            'StateTransitionFcn',      @(x,dt) constturn(x,dt), ... % you provide
            'ProcessNoise',            Qct, ...
            'MeasurementFcn',          @(x)    castLike(Hct6_mat,x)*x, ...
            'MeasurementJacobianFcn',  @(x)    castLike(Hct6_mat,x), ...
            'MeasurementNoise',        R6);

        models    = {ekfCV, ekfCA, ekfCT};
        names     = ["constvel","constacc","constturn"];
        TPM       = [0.986 0.010 0.004;
                     0.010 0.986 0.004;
                     0.015 0.015 0.970];
        initProbs = [0.85 0.10 0.05];
    end

    %% ======= IMM (no unsupported name-value pairs) =======
    imm = trackingIMM( ...
        models, ...                       % filters
        @switchimm, ...                   % your model conversion function
        TPM, ...                          % transition probabilities
        'ModelNames',        cellstr(names), ...
        'ModelProbabilities', cast(initProbs,'like',clsLike), ...
        'MeasurementNoise',  R6);         % optional, OK to include

    % Expect 6-D measurements (positions + velocities)
    setMeasurementSizes(imm, 6, 6);
end
