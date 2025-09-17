%#codegen
classdef tracking_SAFEIMM < matlabshared.tracking.internal.AbstractTrackingFilter ...
        & matlabshared.tracking.internal.fusion.CustomDisplay ...
        & matlabshared.tracking.internal.AbstractJPDAFilter ...
        & fusion.internal.IMMSmoother ...
        & matlabshared.tracking.internal.RetrodictionFilter ...
        & matlabshared.tracking.internal.fusion.AbstractTunableFilter

    %======================== Dependent public get =========================
    properties (Dependent, SetAccess = private, GetAccess = public)
        TrackingFilters
    end

    properties (Access = public)
        ModelConversionFcn
    end

    properties (Dependent)
        TransitionProbabilities
    end

    properties (Dependent, SetAccess = private, GetAccess = public)
        ModelProbabilities
    end

    properties (Dependent, SetAccess = private, GetAccess = public)
        State
        StateCovariance
    end

    properties (Dependent, SetAccess = public)
        MeasurementNoise
    end

    properties (SetAccess = protected, Dependent)
        HasMeasurementWrapping
    end

    properties (SetAccess = protected)
        ModelNames
    end

    %============================= Private data ===========================
    properties (Access = protected)
        pState
        pStateCovariance
        pTrackingFilters
        pTransitionProbabilities         % base Π (kept for compatibility)
        pTPM_Base                        % explicit base Π
        pTPM_Current                     % Π_k (adaptive per-step)
        pMeasurementNoise
        pModelProbabilities
        pNumModels
        pIsLinearKalmanFilter
        pLastRetrodictionDT

        % defaults for dynamic F/Q when using trackingKF('MotionModel','Custom')
        pCVq = single(0.5)
        pCAq = single(0.2)
    end

    properties (Constant, Access = private)
        defaultTransitionProbabilities = single(0.9);
    end

    %==================== BESTCHOICE-IMM =============
    properties (Access = protected)
        % Winner/state memory
        pBC_LastWinner = uint32(1);
        pBC_LastYaw    = single(0);
        pBC_HasVel     = false;
        pBC_LastVel    = zeros(3,1,'single');

        % Heuristic tunables
        pBC_YawRateThresh = single(0.10);
        pBC_AccelThresh   = single(1.50);
        pBC_Bias          = single(0.20);
        pBC_Hysteresis    = single(0.05);
        pBC_StickyDecay   = single(0.30);

        % If true, *try* WTA; in this novel version, a gate decides per-step.
        pBC_UseWinnerTakesAll = true;

        % Constructor NV pickup flag
        pBC_NVParsed = false;
    end

    properties (Access = protected)
        pEPS_CacheValid logical = false
        pEPS_R_chol
        pEPS_trPbar
        pEPS_LastWinner uint32 = uint32(0)
        pEPS_TrPRelTol  single = single(5e-2)   % 5% relative trace tolerance
        pEPS_LastHash   uint32 = uint32(0)

        pSharpMinBF     single = single(2.0);  % log Bayes factor threshold for sharpness
    end

    %==================== Testing===========================================
    properties
        pDBG_Enable logical = false           % turn on/off debug capture
        pDBG_WTA logical = false
        pDBG_EpsGateFired logical = false
        pDBG_EpsBound double = 0
        pDBG_MixWinnerNorm double = 0
        pDBG_JamActive logical = false
        pDBG_CUSUM double = 0
        pDBG_RidgeMult uint16 = 1
        pDBG_UpdatedThisStep logical = false
    end


    %==================== Fast profile controls =====================
    properties (Access = protected)
        % GLR/CUSUM gate for WTA
        pGLR_CUSUM   = single(0);
        pGLR_ThreshA = single(5.0);
        pGLR_Drift   = single(0.10);

        % Safety bound / separation
        pMixTailMax  = single(0.15);
        pMah2Min     = single(9.0);

        % Variable-Structure IMM
        pActiveMask
        pLowCount
        pVS_DeltaOff = single(0.20);
        pVS_DeltaOn  = single(0.40);
        pVS_Toff     = uint16(1);
        pEnsureOneActive = true;

        % Fast profile switches
        pFastMode       = true;
        pUseEntropyCtrl = true;          % Use a cheap sharpness check
        pEntropyTarget  = single(0.25);  % kept for API

        % Inactive model decimation
        pSkipInactive = uint8(3);
        pInactiveTicks

        % Pre-decoded model kind: 1=CV-like, 2=CA-like
        pModelKind
    end

    % Debug/telemetry (optional)
    properties
        pDBG_WTA_Attempts   uint32 = uint32(0);
        pDBG_WTA_Fired      uint32 = uint32(0);
        pDBG_WTA_EpsFail    uint32 = uint32(0);
        pDBG_WTA_SharpFail  uint32 = uint32(0);
        pDBG_WTA_GLRFail    uint32 = uint32(0);
        pDBG_Last_Beps      single = single(0);
        pDBG_Last_EpsLimit  single = single(0);
        pDBG_Last_pmax      single = single(0);
        pDBG_Last_tMass     single = single(0);
        pDBG_WTA_EligibleEPS   uint32 = uint32(0);   % NEW
        pDBG_WTA_EligibleOther uint32 = uint32(0);   % NEW
        pDBG_EpsViolations     uint32 = uint32(0);   % NEW
        pDBG_TailMass  single = single(0);
        pDBG_trPbar    single = single(0);
        pDBG_avgd2     single = single(0);

    end


    %==================== Robustness + Adaptation + Guarantees ========
    properties (Access = protected)
        % --- Robust likelihood selection: 0=None, 1=Student-t, 2=Huber
        pRobustKind      = uint8(1);      % default Student-, for ablation study 0
        pNu              = single(5.0);   % Student-t degrees of freedom
        pNu_Jam          = single(1.0);   % stronger heavy-tail under jam
        pHuberK          = single(3.0);   % Huber threshold (in sqrt(NIS))
        pHuberK_Jam      = single(1.0);   % tighter under jam

        % --- Jamming / clutter detector (CUSUM on NIS)
        pJAM_CUSUM       = single(0);
        pJAM_ThreshA     = single(8.0);
        pJAM_Drift       = single(0.10);
        pJamActive       = true;

        % --- ε-safe tuning 
        pEpsMargin         single = single(0.05);  % require B <= ε - δ
        pEpsHystReq        uint8  = uint8(2);      % frames needed below threshold
        pEpsHystStreak     uint8  = uint8(0);      % internal counter


        % --- Adaptive Π_k controls
        pAdaptAlphaMax   = single(0.7);   % max blending to boosted Π
        pAdaptGLRGain    = single(0.10);  % maps GLR_CUSUM to alpha
        pAdaptEntGain    = single(0.50);  % maps entropy shortfall to alpha
        pAdaptWinnerBias = single(0.10);  % bias out of current winner row
        pBoostToCA       = single(0.15);  % how much to boost towards CA-like
        pBoostToCV       = single(0.05);  % (optionally) towards CV when calm

        % --- ε-safe guarantees
        pEpsSafeWTA      = single(0.50);  % bound for WTA (tail*sqrt(sep)≤ε)
        pEpsInactiveMass = single(0.10);  % never park more than 10% mass
        pTopKMax         = uint8(3);      % fallback to top-K mixing if needed
    end

    %============================== Methods ================================
    methods
        function IMM = tracking_SAFEIMM(varargin)
            % Constructor (preserves external API).

            % Parse smoothing args then remove them for our own parsing
            [smoothArgs, idx1, idx2] = trackingSAFEIMM.parseSmoothingNVPairs(varargin{:});
            vArgsNoSmoothing = {varargin{1:idx1-1}, varargin{idx1+2:idx2-1}, varargin{idx2+2:end}};

            % Parse OOSM args and remove them for filter parsing
            [oosmArgs, idx1] = trackingEKF.parseOOSMNVPairs(vArgsNoSmoothing{:});
            filterArgs = {vArgsNoSmoothing{1:idx1-1}, vArgsNoSmoothing{idx1+2:end}};

            IMM@fusion.internal.IMMSmoother(smoothArgs{:});

            % Pick up heuristic NVs directly (ignore unknown keys)
            if ~IMM.pBC_NVParsed
                bestIdx = matlabshared.tracking.internal.findFirstNVPair(varargin{:});
                if bestIdx <= numel(varargin)
                    nvTail = varargin(bestIdx:end);
                    for k_nv = 1:2:numel(nvTail)
                        key = nvTail{k_nv};
                        if ~(ischar(key) || isstring(key)), continue; end
                        val = nvTail{k_nv+1};
                        switch lower(string(key))
                            case "yawratethresh",     IMM.pBC_YawRateThresh = single(val);
                            case "accelthresh",       IMM.pBC_AccelThresh   = single(val);
                            case "bias",              IMM.pBC_Bias          = single(val);
                            case "hysteresis",        IMM.pBC_Hysteresis    = single(val);
                            case "stickydecay",       IMM.pBC_StickyDecay   = single(val);
                            case "usewinnertakesall", IMM.pBC_UseWinnerTakesAll = logical(val);
                            otherwise
                                % ignore
                        end
                    end
                end
                IMM.pBC_NVParsed = true;
            end

            % Parse constructor inputs (state, filters, TPM, etc.)
            [state, stateCov, trackingFilters, transitionProbabilities, modelConversion, ...
                modelProb, measNoise, modelNames] = parseInputs(IMM, filterArgs{:});

            numModels              = numel(trackingFilters);
            IMM.pNumModels         = numModels;
            IMM.TrackingFilters    = trackingFilters;
            IMM.TransitionProbabilities = transitionProbabilities;   % sets base + current
            IMM.ModelConversionFcn = modelConversion;

            if ~isempty(measNoise), IMM.MeasurementNoise = measNoise; end

            classToUse = class(IMM.pTrackingFilters{1}.State);

            if isempty(modelProb)
                IMM.pModelProbabilities = ones(numModels, 1, 'like', IMM.TransitionProbabilities)/numModels;
            else
                IMM.ModelProbabilities = cast(modelProb, classToUse);
            end

            if isempty(modelNames)
                IMM.ModelNames = defaultModelNames(trackingFilters);
            else
                validateattributes(modelNames, {'cell','string'}, {'numel', numel(trackingFilters)}, 'trackingSAFEIMM', 'ModelNames');
                if isstring(modelNames), modelNames = cellstr(modelNames); end
                IMM.ModelNames = modelNames;
            end

            % Init combined state (convert/mix if not provided)
            m = numel(IMM.pTrackingFilters{1}.State);
            if (isempty(state) && isempty(stateCov))
                [IMM.pState, IMM.pStateCovariance] = IMM.weightedCombine(IMM.pModelProbabilities);
            else
                if (isempty(state) && ~isempty(stateCov))
                    newState    = cast(matlabshared.tracking.internal.expandScalarValue(0, [m,1]), classToUse);
                    newStateCov = cast(matlabshared.tracking.internal.expandScalarValue(stateCov, [m,m]), classToUse);
                elseif (~isempty(state) && isempty(stateCov))
                    newState    = cast(matlabshared.tracking.internal.expandScalarValue(state, [m,1]), classToUse);
                    newStateCov = cast(matlabshared.tracking.internal.expandScalarValue(1, [m,m]), classToUse);
                else
                    newState    = cast(matlabshared.tracking.internal.expandScalarValue(state, [m,1]), classToUse);
                    newStateCov = cast(matlabshared.tracking.internal.expandScalarValue(stateCov, [m,m]), classToUse);
                end
                IMM.State           = newState;
                IMM.StateCovariance = newStateCov;
                initialize(IMM, newState, newStateCov);
            end

            % Setup OOSM memory
            if ~isempty(oosmArgs)
                if oosmArgs{2} > 0
                    IMM.MaxNumOOSMSteps = oosmArgs{2};
                end
            end

            % Initializations
            IMM.pBC_LastWinner = uint32(1);
            [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(IMM.State);

            % Decision-theoretic init
            IMM.pActiveMask    = true(IMM.pNumModels,1);
            IMM.pLowCount      = zeros(IMM.pNumModels,1,'uint16');
            IMM.pInactiveTicks = zeros(IMM.pNumModels,1,'uint8');
            IMM.pGLR_CUSUM     = single(0);
            IMM.pJAM_CUSUM     = single(0);
            IMM.pJamActive     = true;

            % Pre-decode model kinds (avoid strings in hot path)
            IMM.pModelKind = zeros(IMM.pNumModels,1,'uint8');
            for k = 1:IMM.pNumModels
                mk = uint8(2); % default CA-like
                if ~isempty(IMM.ModelNames)
                    nm = lower(string(IMM.ModelNames{k}));
                    if startsWith(nm,"cv") || contains(nm,"constvel")
                        mk = uint8(1);
                    end
                end
                IMM.pModelKind(k) = mk;
            end

            % Initialize adaptive Π_k
            IMM.pTPM_Current = IMM.pTPM_Base;
        end

        %------------------------------ Setters ---------------------------
        function set.State(IMM, value)
            validateattributes(value, {'double','single'}, {'real','finite','nonsparse','vector'}, 'trackingSAFEIMM', 'State');
            if coder.internal.is_defined(IMM.pTrackingFilters{1}.State)
                m          = numel(IMM.pTrackingFilters{1}.State);
                classToUse = class(IMM.pTrackingFilters{1}.State);
                value      = cast(value, classToUse);
                coder.internal.assert(numel(value) == m, 'fusion:trackingSAFEIMM:invalidStateDims', 'State', m);
            end
            IMM.pState = value(:);
        end

        function value = get.State(IMM),          value = IMM.pState; end

        function set.StateCovariance(IMM, value)
            validateattributes(value, {'double','single'}, {'real','finite','nonsparse','2d','nonempty','square'}, 'trackingSAFEIMM', 'StateCovariance');
            if coder.internal.is_defined(IMM.pTrackingFilters{1}.State)
                m          = numel(IMM.pTrackingFilters{1}.State);
                classToUse = class(IMM.pTrackingFilters{1}.State);
                value      = cast(value, classToUse);
                coder.internal.assert(size(value,1) == m, 'fusion:trackingSAFEIMM:invalidStateCovDims', 'StateCovariance', m);
            end
            matlabshared.tracking.internal.isSymmetricPositiveSemiDefinite('StateCovariance', value);
            IMM.pStateCovariance = value;
        end

        function value = get.StateCovariance(IMM), value = IMM.pStateCovariance; end

        function set.TrackingFilters(IMM, value)
            validateattributes(value, {'cell'}, {'nonempty'}, 'trackingSAFEIMM', 'TrackingFilters');
            coder.internal.errorIf(IMM.pNumModels <= 1, 'fusion:trackingSAFEIMM:invalidMultiModel', 'TrackingFilters');
            IMM.pTrackingFilters      = cell(IMM.pNumModels, 1);
            isLinearFilter            = false(1, IMM.pNumModels);
            for k = coder.unroll(1:IMM.pNumModels)
                isMMCompatible = fusion.internal.isMultiModelCompatible(value{k});
                coder.internal.errorIf(~isMMCompatible(2), 'fusion:trackingSAFEIMM:invalidFilterType', class(value{k}));
                cond = ~isequal(class(value{1}.State), class(value{k}.State));
                coder.internal.errorIf(cond, 'fusion:trackingSAFEIMM:invalidStateType', 'TrackingFilters', class(value{k}.State));
                IMM.pTrackingFilters{k} = clone(value{k});
                isLinearFilter(k)       = matlabshared.tracking.internal.isLinearKalmanFilter(IMM.pTrackingFilters{k});
            end
            IMM.pIsLinearKalmanFilter = isLinearFilter;
        end

        function val = get.TrackingFilters(IMM), val = IMM.pTrackingFilters; end

        function set.TransitionProbabilities(IMM, value)
            validateattributes(value, {'double','single'}, {'real','finite','nonsparse','nonempty','nonnegative','<=',1}, 'trackingSAFEIMM', 'TransitionProbabilities');
            classToUse = class(IMM.pTrackingFilters{1}.State);
            numModels  = IMM.pNumModels;
            if isscalar(value)
                offDiag  = cast((1 - value)/(numModels - 1), classToUse);
                transProb = offDiag*ones(numModels, classToUse);
                idx       = 1:(numModels + 1):numModels*numModels;
                transProb(idx) = value;
            elseif isvector(value)
                coder.internal.errorIf(~isequal(numel(value), IMM.pNumModels), 'fusion:trackingSAFEIMM:invalidNumProbVector', 'TransitionProbabilities', IMM.pNumModels);
                transProb = zeros(numModels, classToUse);
                for j = 1:numModels
                    transProb(j,:) = (1 - value(j))/(numModels - 1);
                    transProb(j,j) = value(j);
                end
            else
                validateattributes(value, {'double','single'}, {'2d','square'}, 'trackingSAFEIMM', 'TransitionProbabilities');
                coder.internal.errorIf(~isequal(size(value,1), IMM.pNumModels), 'fusion:trackingSAFEIMM:invalidNumProbMatrix', 'TransitionProbabilities', IMM.pNumModels);
                cond = ~(sum(abs(sum(value,2) - ones(IMM.pNumModels,1,classToUse))) < cast(1e-5, classToUse));
                coder.internal.errorIf(cond, 'fusion:trackingSAFEIMM:invalidTransProbMatrix', 'TransitionProbabilities');
                transProb = value;
            end
            IMM.pTransitionProbabilities = cast(transProb, classToUse);
            IMM.pTPM_Base = IMM.pTransitionProbabilities;  % keep base
            IMM.pTPM_Current = IMM.pTPM_Base;             % initialize current
        end

        function val = get.TransitionProbabilities(IMM), val = IMM.pTransitionProbabilities; end

        function set.ModelConversionFcn(IMM, value)
            validateattributes(value, {'function_handle'}, {'nonempty'}, 'trackingSAFEIMM', 'ModelconversionFcn');
            IMM.ModelConversionFcn = value;
        end
        function val = get.ModelConversionFcn(IMM), val = IMM.ModelConversionFcn; end

        function set.ModelProbabilities(IMM, value)
            validateattributes(value, {'double','single'}, {'real','finite','nonsparse','nonempty','vector','numel',IMM.pNumModels}, 'trackingSAFEIMM', 'ModelProbabilities');
            value                 = value/sum(value);
            IMM.pModelProbabilities = value;
        end
        function val = get.ModelProbabilities(IMM), val = IMM.pModelProbabilities; end

        function set.pModelProbabilities(IMM, val)
            if coder.internal.is_defined(IMM.pModelProbabilities)
                IMM.pModelProbabilities(:) = cast(val, 'like', IMM.pModelProbabilities);
            else
                IMM.pModelProbabilities = val;
            end
        end

        function set.MeasurementNoise(IMM, value)
            validateattributes(value, {'single','double'}, {'real','finite','nonsparse','2d','nonempty','square'}, 'trackingSAFEIMM', 'MeasurementNoise');
            classToUse = class(IMM.pTrackingFilters{1}.State);
            matlabshared.tracking.internal.isSymmetricPositiveSemiDefinite('MeasurementNoise', value);
            for i = coder.unroll(1:IMM.pNumModels)
                IMM.pTrackingFilters{i}.MeasurementNoise = cast(value, classToUse);
            end
        end
        function value = get.MeasurementNoise(IMM), value = IMM.pTrackingFilters{1}.MeasurementNoise; end

        function value = get.HasMeasurementWrapping(IMM)
            value = zeros(1, IMM.pNumModels, 'logical');
            for i = coder.unroll(1:IMM.pNumModels)
                if isprop(IMM.pTrackingFilters{i}, 'HasMeasurementWrapping')
                    value(i) = IMM.pTrackingFilters{i}.HasMeasurementWrapping;
                end
            end
        end

        %------------------------------ initialize -------------------------
        function initialize(IMM, state, stateCov, varargin)
            coder.internal.errorIf(numel(IMM) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'initialize');
            stateSize = numel(IMM.pTrackingFilters{1}.State);
            index     = (1:stateSize);
            s         = state(index');
            sc        = stateCov(1:stateSize, 1:stateSize);
            givenModel = 1;

            if isa(IMM.TrackingFilters{givenModel}, 'trackingPF')
                initialize(IMM.TrackingFilters{givenModel}, IMM.TrackingFilters{givenModel}.NumParticles, s, sc);
            else
                IMM.TrackingFilters{givenModel}.State           = s;
                IMM.TrackingFilters{givenModel}.StateCovariance = sc;
            end

            for mdlOut = coder.unroll(2:IMM.pNumModels)
                Xzeros = zeros(size(IMM.TrackingFilters{mdlOut}.State), 'like', IMM.TrackingFilters{mdlOut}.State);
                Pzeros = zeros(size(IMM.TrackingFilters{mdlOut}.StateCovariance), 'like', IMM.TrackingFilters{mdlOut}.StateCovariance);
                Xk = IMM.ModelConversionFcn(IMM.ModelNames{givenModel}, s,  IMM.ModelNames{mdlOut}, Xzeros);
                if nargin(IMM.ModelConversionFcn) == 5
                    Pk = IMM.ModelConversionFcn(IMM.ModelNames{givenModel}, sc, IMM.ModelNames{mdlOut}, Pzeros, s);
                else
                    Pk = IMM.ModelConversionFcn(IMM.ModelNames{givenModel}, sc, IMM.ModelNames{mdlOut}, Pzeros);
                end
                if isa(IMM.TrackingFilters{mdlOut}, 'trackingPF')
                    initialize(IMM.TrackingFilters{mdlOut}, IMM.TrackingFilters{mdlOut}.NumParticles, Xk, Pk);
                else
                    IMM.TrackingFilters{mdlOut}.State           = Xk;
                    IMM.TrackingFilters{mdlOut}.StateCovariance = Pk;
                end
            end

            [IMM.pState, IMM.pStateCovariance] = IMM.weightedCombine(IMM.pModelProbabilities);

            coder.internal.assert(mod(numel(varargin),2) == 0, 'MATLAB:system:invalidPVPairs');
            for i = 1:2:numel(varargin)
                param = validatestring(varargin{i}, {'ModelProbabilities','TransitionProbabilities','MeasurementNoise'}, mfilename);
                IMM.(param) = varargin{i+1};
            end
        end

        %------------------------------ predict ----------------------------
        function [x_pred, P_pred] = predict(IMM, varargin)
            coder.internal.errorIf(numel(IMM) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'predict');

            if nargin < 2
                dt = ones(1, 'like', IMM.pState);
            else
                validateattributes(varargin{1}, {'double','single'}, {'real','finite','nonsparse','nonnegative','scalar'}, 'predict', 'dt');
                dt = varargin{1};
            end

            % Adaptive Π_k update from last correction info (GLR, entropy)
            IMM.updateAdaptiveTPM();

            % Priors: c <- Π_k' * c (dynamic TPM)
            IMM.advancePriors();  % Uses pTPM_Current

            % Predict (decimate inactive)
            for kMdl = coder.unroll(1:IMM.pNumModels)
                if ~IMM.pActiveMask(kMdl)
                    if IMM.pInactiveTicks(kMdl) < IMM.pSkipInactive
                        IMM.pInactiveTicks(kMdl) = IMM.pInactiveTicks(kMdl) + 1;
                        continue;
                    else
                        IMM.pInactiveTicks(kMdl) = uint8(0);
                    end
                end

                fk = IMM.pTrackingFilters{kMdl};
                if isa(fk,'trackingKF') && isprop(fk,'MotionModel') && strcmpi(fk.MotionModel,'Custom')
                    if IMM.pModelKind(kMdl) == 1 % CV-like
                        fk.StateTransitionModel = IMM.bcMakeFcv(dt, fk.State);
                        fk.ProcessNoise         = IMM.bcMakeQcv(dt, fk.State);
                    else                        % CA-like
                        fk.StateTransitionModel = IMM.bcMakeFca(dt, fk.State);
                        fk.ProcessNoise         = IMM.bcMakeQca(dt, fk.State);
                    end
                    predict(fk);
                else
                    predict(fk, dt);
                end
            end

            % Output prediction from last winner for continuity
            [x_pred, P_pred] = IMM.readoutAfterPredict();
            IMM.LastStepTime = cast(dt, 'like', IMM.LastStepTime);
            updatePredictionData(IMM);
            updateRetrodictionStateAfterPrediction(IMM, dt);
        end


        %------------------------------ correct ----------------------------
        function [x_corr, P_corr, mdl_probs] = correct(IMM, z, varargin)
            narginchk(2, inf);
            coder.internal.errorIf(numel(IMM) > 1, ...
                'fusion:trackingSAFEIMM:NonScalarFilter', ...
                'trackingSAFEIMM', 'correct');

            % ---- guard: missed or invalid measurement -> coast (no state change) ----
            if isempty(z) || ~all(isfinite(z(:)))
                mdl_probs = IMM.pModelProbabilities(:);
                [x_corr, P_corr] = IMM.readoutAfterPredict();
                updateCorrectionData(IMM);
                updateHistoryAfterCorrection(IMM);
                x_corr = double(x_corr); P_corr = double(P_corr);
                return
            end


            % Keep all numerics in the same class as model probabilities
            likeP       = IMM.pModelProbabilities;
            Mdim        = max(size(z,1), size(z,2));                   % robust for row/col
            likelihoods = zeros(IMM.pNumModels, 1, 'like', likeP);
            nis_vec     = zeros(IMM.pNumModels, 1, 'like', likeP);

            % Flag: did any model actually apply a correction this frame?
            updatesApplied = false;

            % --- Outlier gate: protect state from extreme residuals under jam ---
            % --- Chi-square(0.999) gate once (Wilson–Hilferty) ---
            m    = cast(max(1, Mdim), 'like', likeP);
            z999 = cast(3.09, 'like', likeP);                    % ≈ Φ^{-1}(0.999)
            term = (cast(1,'like',likeP) - 2./(9*m) + z999.*sqrt(2./(9*m)));
            d2GateNom = m .* term .* term .* term;               % ~chi2inv(0.999, m)


            % Under jam: relax the gate a bit (codegen-safe scalars)
            d2Gate = d2GateNom;
            if IMM.pJamActive
                d2Gate = cast(2.0,'like',d2GateNom) * d2GateNom;   % relax under jam
            end

            % Per-model likelihood & (possibly) correction
            for kMdl = coder.unroll(1:IMM.pNumModels)
                isPrevWinner = (uint32(kMdl) == IMM.pBC_LastWinner);
                if ~IMM.pActiveMask(kMdl) && ~isPrevWinner
                    likelihoods(kMdl) = eps('like', likeP);
                    nis_vec(kMdl)     = cast(0, 'like', likeP);
                    continue;
                end


                fk = IMM.pTrackingFilters{kMdl};

                % Compute NIS (codegen-safe) and keep class
                d2   = IMM.safeDistance(fk, z, varargin{:});
                d2   = max(d2, 0);
                d2c  = cast(d2, 'like', likeP);
                nis_vec(kMdl) = d2c;

                % --- Likelihood family (one per step to avoid scale mismatch) ---

                lg    = likelihood(fk, z, varargin{:});                 % normalized, scale-aware
                lg    = max(lg, realmin('like', lg));
                if IMM.pJamActive
                    lr    = cast(IMM.robustLikelihoodFromNIS(d2, Mdim), 'like', lg);
                    kg    = exp(-0.5 * cast(d2, 'like', lg));        % no double promotion
                    tiny  = realmin('like', lg);
                    ratio = lr / max(kg, tiny);
                    % clamp ratio to avoid overflow/underflow in extreme jams
                    ratio = min(max(ratio, cast(1e-12,'like',lg)), cast(1e12,'like',lg));
                    likelihoods(kMdl) = max(lg, tiny) * ratio;
                else
                    likelihoods(kMdl) = lg;
                end      

                doUpdate     = (d2c <= d2Gate) || isPrevWinner;         
                if doUpdate
                    % ---- State update (R-inflated only under jam) ----
                    if IMM.pJamActive && isprop(fk, 'MeasurementNoise')
                        R_saved = fk.MeasurementNoise;
                        scaleR  = cast(1.0, 'like', R_saved);            % 2–3× is typical
                        fk.MeasurementNoise = scaleR * R_saved;

                        correct(fk, z, varargin{:});
                        % debug command
                        IMM.pDBG_UpdatedThisStep = true;

                        fk.MeasurementNoise = R_saved;                   % restore
                    else
                        correct(fk, z, varargin{:});
                        % debug command
                        IMM.pDBG_UpdatedThisStep = true;
                    end
                    updatesApplied = true;
                end


            end

            if ~updatesApplied
                % Force a correction on the previous winner so the track never stalls
                kForce  = double(IMM.pBC_LastWinner);
                fkForce = IMM.pTrackingFilters{kForce};
                correct(fkForce, z, varargin{:});
            end


            % --- Posterior over models (log-domain mixing; prevents underflow) ---
            tiny = realmin('like', likeP);
            pr   = max(IMM.pModelProbabilities(:), tiny);
            lk   = max(likelihoods(:),           tiny);

            L  = log(lk) + log(pr);                  % same class as likeP
            Lm = max(L);
            w  = exp(L - Lm);
            s  = sum(w);
            if s > 0
                mdl_probs = w / s;
            else
                mdl_probs = ones(IMM.pNumModels, 1, 'like', w) / IMM.pNumModels;
            end

            % --- Jamming CUSUM: normalize by expected NIS (= Mdim) ---
            avgNIS = sum(mdl_probs .* nis_vec);
            IMM.updateJammingCUSUM(avgNIS, Mdim);

            % Winner & gates
            [pmax, widx, p2] = IMM.bestTwoWithValues(mdl_probs);
            tinyLog = tiny;
            if p2 == 0
                logBF = cast(inf, 'like', likeP);
            else
                logBF = log(max(pmax, tinyLog)) - log(max(p2, tinyLog));
            end
            IMM.pGLR_CUSUM = max(cast(0,'like',likeP), IMM.pGLR_CUSUM + cast(logBF,'like',IMM.pGLR_CUSUM) - IMM.pGLR_Drift);

            % --- softer GLR & sharpness gates (local scalars; no property writes) ---
            thA = cast(0.6,'like',likeP) * IMM.pGLR_ThreshA;   % 40% softer than default
            passGLR = (IMM.pGLR_CUSUM >= thA);

            % simple sharpness on probabilities (avoids brittle Bayes-factor only)
            passSharp = (pmax >= cast(0.60,'like',pmax));      % ~60% winner mass

            % --- compute ε on the SAME Top-K mixture you’d output (use K<=2 to keep tail small) ---
            Kcap    = uint8( min( double(IMM.pTopKMax), 2 ) );
            keepIdx = IMM.selectTopK(mdl_probs, Kcap);
            weights_eps = zeros(IMM.pNumModels,1,'like',mdl_probs);
            for jj = 1:numel(keepIdx)
                k = keepIdx(jj);
                if k>=1 && k<=IMM.pNumModels
                    weights_eps(k) = mdl_probs(k);
                end
            end

            % B_eps    = IMM.epsSafeBoundWTA(widx, weights_eps);
            % %disp(B_eps)
            % epsLimit = cast(2.0,'like',B_eps);                 % start with 2.0; you can lower later
            % epsOK    = (B_eps <= epsLimit);
            % 
            % % --- final decision: ε must pass AND (posterior sharp OR GLR high) ---
            % doWTA = IMM.pFastMode && IMM.pBC_UseWinnerTakesAll && epsOK && (passSharp || passGLR);

            B_eps    = IMM.epsSafeBoundWTA(widx, weights_eps);
            epsLimit = cast(2.0,'like',B_eps);     % keep your existing limit
            epsOKraw = (B_eps <= epsLimit);

            % Eligibility counters
            if epsOKraw
                IMM.pDBG_WTA_EligibleEPS = IMM.pDBG_WTA_EligibleEPS + uint32(1);
            end
            if (passSharp || passGLR)
                IMM.pDBG_WTA_EligibleOther = IMM.pDBG_WTA_EligibleOther + uint32(1);
            end

            % ε-margin and hysteresis
            epsOKmargin = (B_eps <= (epsLimit - IMM.pEpsMargin));
            if epsOKmargin
                IMM.pEpsHystStreak = uint8(min(255, IMM.pEpsHystStreak + 1));
            else
                IMM.pEpsHystStreak = uint8(0);
            end
            epsOK = epsOKmargin && (IMM.pEpsHystStreak >= IMM.pEpsHystReq);

            % Final decision: ε must pass + (sharp OR GLR)
            doWTA = IMM.pFastMode && IMM.pBC_UseWinnerTakesAll && epsOK && (passSharp || passGLR);
              
            % disable SAFE, for ablation study 
            % doWTA = IMM.pBC_UseWinnerTakesAll;

            

            IMM.pDBG_WTA_Attempts = IMM.pDBG_WTA_Attempts + uint32(1);
            IMM.pDBG_Last_Beps     = single(B_eps);
            IMM.pDBG_Last_EpsLimit = single(epsLimit);
            IMM.pDBG_Last_pmax     = single(pmax);
            IMM.pDBG_Last_tMass    = single(1 - pmax);  % crude tail proxy
            IMM.pDBG_TailMass = single(1 - pmax);
            IMM.pDBG_EpsGateFired = logical(epsOK);   % ε condition (with margin+hysteresis)


            if doWTA
                IMM.pDBG_WTA_Fired = IMM.pDBG_WTA_Fired + uint32(1);
                IMM.pDBG_WTA = true;
            else
                IMM.pDBG_WTA = false;
                if ~epsOK,     IMM.pDBG_WTA_EpsFail   = IMM.pDBG_WTA_EpsFail   + uint32(1); end
                if ~passSharp, IMM.pDBG_WTA_SharpFail = IMM.pDBG_WTA_SharpFail + uint32(1); end
                if ~passGLR,   IMM.pDBG_WTA_GLRFail   = IMM.pDBG_WTA_GLRFail   + uint32(1); end
            end


            % debug and audit
            if IMM.pDBG_Enable
                % Winner state in winner's space
                [xw_audit, ~] = IMM.readFromModel(widx);
                dw = size(xw_audit,1);
                M  = IMM.pNumModels;

                % Top-K indices (may be variable-length)
                keepIdxAudit = IMM.selectTopK(mdl_probs, IMM.pTopKMax);

                % Build a fixed-size selection mask for codegen
                sel = false(M,1);
                for j = 1:numel(keepIdxAudit)
                    idx = keepIdxAudit(j);
                    if idx >= 1 && idx <= M
                        sel(idx) = true;
                    end
                end

                % Sum of selected weights (normalize safely)
                wsel_sum = 0;
                for i = 1:M
                    if sel(i)
                        wsel_sum = wsel_sum + mdl_probs(i);
                    end
                end
                if wsel_sum <= 0
                    wsel_sum = 1;  % avoid divide-by-zero; yields zeros mix below
                end

                % Mixture mean expressed in the WINNER's space
                xmix_w = zeros(dw,1,'like',xw_audit);
                for i = 1:M
                    if sel(i)
                        fi = IMM.pTrackingFilters{i};
                        [mu_i_w, ~] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, widx);
                        % Guard: only accumulate if dimensions match winner space
                        if isequal(size(mu_i_w), size(xw_audit))
                            xmix_w = xmix_w + (mdl_probs(i)/wsel_sum) * mu_i_w;
                        end
                    end
                end

                % Euclidean distance in a numeric-stable way (no try/catch)
                diff = double(xmix_w) - double(xw_audit);
                IMM.pDBG_MixWinnerNorm = sqrt(sum(diff.^2));
            end



            % --- Safer model activation during jam (temporary raise of OFF bar) ---
            % --- Jam-time activation (preserve your result) with skip-proof + ε cap ---
            if IMM.pJamActive
                % 1) Your original “raise OFF bar” update
                oldOff = IMM.pVS_DeltaOff;
                IMM.pVS_DeltaOff = max(oldOff, cast(0.25,'like',oldOff));
                IMM.updateActiveMask(mdl_probs, widx);
                IMM.pVS_DeltaOff = oldOff;

                % 2) Always keep winner active
                IMM.pActiveMask(widx) = true;

                % 3) ε-safe cap on inactive mass (reactivate strongest until under cap)
                inactMass = sum(mdl_probs(~IMM.pActiveMask));
                if inactMass > IMM.pEpsInactiveMass
                    idxInact = find(~IMM.pActiveMask);
                    % sort inactive by descending probability (codegen-safe)
                    for a = 1:numel(idxInact)-1
                        for b = a+1:numel(idxInact)
                            if mdl_probs(idxInact(b)) > mdl_probs(idxInact(a))
                                t = idxInact(a); idxInact(a) = idxInact(b); idxInact(b) = t;
                            end
                        end
                    end
                    k = 1;
                    while (inactMass > IMM.pEpsInactiveMass) && (k <= numel(idxInact))
                        IMM.pActiveMask(idxInact(k)) = true;
                        inactMass = sum(mdl_probs(~IMM.pActiveMask));
                        k = k + 1;
                    end
                end

                % 4) NO-SKIP GUARANTEE for next predict:
                %    For every INACTIVE model, force tick to pSkipInactive so the
                %    predict() gate takes the "do predict" branch once, not 'continue'.
                tickMax = IMM.pSkipInactive;     % uint8
                for i = 1:IMM.pNumModels
                    if IMM.pActiveMask(i)
                        IMM.pInactiveTicks(i) = uint8(0);
                        IMM.pLowCount(i)      = uint16(0);  % avoid hidden carry-over
                    else
                        IMM.pInactiveTicks(i) = tickMax;    % will be predicted next step
                    end
                end
                % --- debug reset (cheap; only used if pDBG_Enable) ---
                IMM.pDBG_UpdatedThisStep = false;
            else
                IMM.updateActiveMask(mdl_probs, widx);
                % Keep winner active
                IMM.pActiveMask(widx) = true;
                % ε-safe cap on inactive mass (reactivate strongest until under cap)
                inactMass = sum(mdl_probs(~IMM.pActiveMask));
                if inactMass > IMM.pEpsInactiveMass
                    idxInact = find(~IMM.pActiveMask);
                    for a = 1:numel(idxInact)-1
                        for b = a+1:numel(idxInact)
                            if mdl_probs(idxInact(b)) > mdl_probs(idxInact(a))
                                t = idxInact(a); idxInact(a) = idxInact(b); idxInact(b) = t;
                            end
                        end
                    end
                    k = 1;
                    while (inactMass > IMM.pEpsInactiveMass) && (k <= numel(idxInact))
                        IMM.pActiveMask(idxInact(k)) = true;
                        inactMass = sum(mdl_probs(~IMM.pActiveMask));
                        k = k + 1;
                    end
                end
            end


            % Debug
            % --- ε-bound audit: only when WTA fires and debug enabled
            if IMM.pDBG_Enable && doWTA
                % slack to avoid FP jitter
                slack = 1e-6;
                if (double(IMM.pDBG_MixWinnerNorm) - double(IMM.pDBG_Last_EpsLimit)) > slack
                    IMM.pDBG_EpsViolations = IMM.pDBG_EpsViolations + uint32(1);
                end
                % expose bound used for this step to logs
                IMM.pDBG_EpsBound = double(IMM.pDBG_Last_EpsLimit);
            end



            % --- Output state (soft WTA to avoid total collapse) ---
            if doWTA
                [IMM.pState, IMM.pStateCovariance] = IMM.readFromModel(widx);
                % Soft one-hot: keep tiny mass on rivals
                IMM.pModelProbabilities = mdl_probs;   % keep soft mass; WTA is readout-only     
                IMM.pBC_LastWinner = uint32(widx);
                [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(IMM.pState);
                x_corr = IMM.pState;
                P_corr = IMM.pStateCovariance;
                
                % This is for checking WTA, HARD WTA
                % [IMM.pState, IMM.pStateCovariance] = IMM.readFromModel(widx);
                % % HARD WTA: propagate a one-hot posterior so priors diverge next step
                % IMM.pModelProbabilities = IMM.oneHot(widx);
                % IMM.pBC_LastWinner = uint32(widx);
                % [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(IMM.pState);
                % x_corr = IMM.pState;
                % P_corr = IMM.pStateCovariance;
            else
                keepIdx = IMM.selectTopK(mdl_probs, IMM.pTopKMax);
                [xmix, Pmix] = IMM.weightedCombine(mdl_probs, keepIdx);
                IMM.pState = xmix; 
                IMM.pStateCovariance = Pmix;
                IMM.pModelProbabilities = mdl_probs;
                IMM.pBC_LastWinner = uint32(widx);
                [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(xmix);
                x_corr = xmix; 
                P_corr = Pmix;
            end
            
            % Debug
            if doWTA
                IMM.pDBG_WTA_Fired = IMM.pDBG_WTA_Fired + uint32(1);
                IMM.pDBG_WTA = true;
            else
                IMM.pDBG_WTA = false;
                if ~epsOK
                    IMM.pDBG_WTA_EpsFail = IMM.pDBG_WTA_EpsFail + uint32(1);
                elseif ~passSharp
                    IMM.pDBG_WTA_SharpFail = IMM.pDBG_WTA_SharpFail + uint32(1);
                elseif ~passGLR
                    IMM.pDBG_WTA_GLRFail = IMM.pDBG_WTA_GLRFail + uint32(1);
                end
            end


            updateCorrectionData(IMM);
            updateHistoryAfterCorrection(IMM);
        end



        %------------------------------ correctjpda ------------------------
        function [x_corr, P_corr, mdl_probs] = correctjpda(IMM, z, jpda, varargin)
            narginchk(3, inf)
            validateCorrectJPDAInputs(IMM, z, jpda, 'correctjpda', varargin{:});

            M = size(z,2);
            mDim = size(z,1);
            likelihoods = zeros(IMM.pNumModels, M, 'like', IMM.pModelProbabilities);
            nis_mat     = zeros(IMM.pNumModels, M, 'like', IMM.pModelProbabilities);

            for kMdl = coder.unroll(1:IMM.pNumModels)
                if ~IMM.pActiveMask(kMdl)
                    likelihoods(kMdl,:) = eps('like', IMM.pModelProbabilities);
                    nis_mat(kMdl,:)     = single(0);
                    continue;
                end
                filterk = IMM.pTrackingFilters{kMdl};
                for meas = 1:M
                    d2 = IMM.safeDistance(filterk, z(:,meas), varargin{:});
                    d2 = max(d2, 0);
                    nis_mat(kMdl,meas) = d2;
                    likelihoods(kMdl,meas) = IMM.robustLikelihoodFromNIS(d2, mDim);
                end
                if IMM.pIsLinearKalmanFilter(kMdl)
                    correctjpda(filterk, z, jpda);
                else
                    correctjpda(filterk, z, jpda, varargin{:});
                end
            end

            mdl_probs = updateModelProbJPDA(IMM, likelihoods, jpda, M);

            % Update jamming CUSUM with JPDA-averaged NIS (normalize by mDim)
            avgNIS = sum(mdl_probs .* (likelihoodsToAvgNIS(nis_mat, jpda)));
            IMM.updateJammingCUSUM(avgNIS, mDim);

            [pmax, widx, p2] = IMM.bestTwoWithValues(mdl_probs);
            if p2==0, logBF = inf; else, logBF = log(max(pmax, realmin)) - log(max(p2, realmin)); end
            IMM.pGLR_CUSUM = max(single(0), IMM.pGLR_CUSUM + single(logBF) - IMM.pGLR_Drift);

            passGLR   = IMM.pGLR_CUSUM >= IMM.pGLR_ThreshA;
            passSharp = IMM.posteriorIsSharp(pmax, p2);

            % --- ε-safe WTA bound ---
            B_eps = IMM.epsSafeBoundWTA(widx, mdl_probs);
            epsSafeOK = (B_eps <= IMM.pEpsSafeWTA);

            doWTA = IMM.pFastMode && IMM.pBC_UseWinnerTakesAll && passSharp && passGLR && epsSafeOK;

            IMM.updateActiveMaskSafe(mdl_probs, widx);

            if doWTA
                [x_corr, P_corr] = IMM.readFromModel(widx);
                IMM.pModelProbabilities = IMM.oneHot(widx);
                IMM.pBC_LastWinner = uint32(widx);
                [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(x_corr);
            else
                keepIdx = IMM.selectTopK(mdl_probs, IMM.pTopKMax);
                [xmix, Pmix] = IMM.weightedCombine(mdl_probs, keepIdx);
                IMM.pState = xmix; IMM.pStateCovariance = Pmix;
                IMM.pModelProbabilities = mdl_probs;
                IMM.pBC_LastWinner = uint32(widx);
                [IMM.pBC_LastVel, IMM.pBC_LastYaw, IMM.pBC_HasVel] = IMM.bcExtractVelYaw(xmix);
                x_corr = xmix; P_corr = Pmix;
            end

            updateCorrectionData(IMM);
            updateHistoryAfterCorrection(IMM);

            function v = likelihoodsToAvgNIS(nis_mat_in, jpda_in)
                % Mix NIS via JPDA association probabilities (simple average)
                beta = jpda_in(:);
                Mloc = size(nis_mat_in,2);
                if numel(beta) >= Mloc
                    v = sum(nis_mat_in .* reshape(beta(1:Mloc).', 1, Mloc), 2);
                else
                    v = mean(nis_mat_in,2);
                end
            end
        end


        %------------------------------ distance --------------------------
        function d = distance(IMM, zMatrix, measParams)
            if nargin < 3, measParams = {}; end
            coder.internal.errorIf((numel(IMM) > 1), 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'distance');
            modelProb = IMM.pModelProbabilities;
            filterk   = IMM.pTrackingFilters{1};
            wk        = modelProb(1);
            if IMM.pIsLinearKalmanFilter(1)
                d = wk*distance(filterk, zMatrix);
            else
                d = wk*distance(filterk, zMatrix, measParams);
            end
            for kMdl = coder.unroll(2:IMM.pNumModels)
                filterk = IMM.pTrackingFilters{kMdl};
                wk      = modelProb(kMdl);
                if IMM.pIsLinearKalmanFilter(kMdl)
                    d = d + wk*distance(filterk, zMatrix);
                else
                    d = d + wk*distance(filterk, zMatrix, measParams);
                end
            end
        end

        %------------------------------ likelihood ------------------------
        function lhood = likelihood(IMM, z, varargin)
            coder.internal.errorIf(numel(IMM) > 1, ...
                'fusion:trackingSAFEIMM:NonScalarFilter', ...
                'trackingSAFEIMM','likelihood');
            validateattributes(z, {'numeric'}, {'real','finite','nonsparse','vector'}, 'likelihood', 'z');

            indLhood = zeros(IMM.pNumModels,1,'like',IMM.State);
            coder.unroll();
            for kMdl = 1:IMM.pNumModels
                if ~IMM.pActiveMask(kMdl)
                    indLhood(kMdl) = eps('like',IMM.State);
                    continue;
                end
                fk = IMM.pTrackingFilters{kMdl};
                % no redundant branch; always forward varargin
                indLhood(kMdl) = likelihood(fk, z, varargin{:});
            end
            cbar  = IMM.pModelProbabilities;
            lhood = sum(cbar(:) .* indLhood);
        end


        %------------------------------ clone -----------------------------
        function newIMM = clone(IMM)
            coder.inline('never');
            coder.internal.errorIf(numel(IMM) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'clone');
            obj        = str2func(coder.const(class(IMM)));
            newFilters = cell(IMM.pNumModels,1);
            for i = coder.unroll(1:IMM.pNumModels)
                newFilters{i} = clone(IMM.pTrackingFilters{i});
            end
            newIMM = obj('TrackingFilters', newFilters, ...
                         'ModelConversionFcn', IMM.ModelConversionFcn, ...
                         'TransitionProbabilities', IMM.TransitionProbabilities, ...
                         'ModelNames', IMM.ModelNames);
            ppProperties = coder.const({'pState','pStateCovariance', ...
                                        'pMeasurementNoise','pModelProbabilities', ...
                                        'pNumModels','pLastRetrodictionDT'});
            for kk = coder.unroll(1:numel(ppProperties))
                if coder.internal.is_defined(IMM.(ppProperties{kk}))
                    newIMM.(ppProperties{kk}) = IMM.(ppProperties{kk});
                end
            end
            copySmootherProperties(IMM, newIMM);
            cloneRetroFilter(newIMM, IMM);

            % Copy tunables
            newIMM.pBC_LastWinner        = IMM.pBC_LastWinner;
            newIMM.pBC_LastYaw           = IMM.pBC_LastYaw;
            newIMM.pBC_HasVel            = IMM.pBC_HasVel;
            newIMM.pBC_LastVel           = IMM.pBC_LastVel;
            newIMM.pBC_YawRateThresh     = IMM.pBC_YawRateThresh;
            newIMM.pBC_AccelThresh       = IMM.pBC_AccelThresh;
            newIMM.pBC_Bias              = IMM.pBC_Bias;
            newIMM.pBC_Hysteresis        = IMM.pBC_Hysteresis;
            newIMM.pBC_StickyDecay       = IMM.pBC_StickyDecay;
            newIMM.pBC_UseWinnerTakesAll = IMM.pBC_UseWinnerTakesAll;
            newIMM.pBC_NVParsed          = IMM.pBC_NVParsed;

            newIMM.pGLR_CUSUM   = IMM.pGLR_CUSUM;
            newIMM.pGLR_ThreshA = IMM.pGLR_ThreshA;
            newIMM.pGLR_Drift   = IMM.pGLR_Drift;
            newIMM.pMixTailMax  = IMM.pMixTailMax;
            newIMM.pMah2Min     = IMM.pMah2Min;
            newIMM.pActiveMask  = IMM.pActiveMask;
            newIMM.pLowCount    = IMM.pLowCount;
            newIMM.pVS_DeltaOff = IMM.pVS_DeltaOff;
            newIMM.pVS_DeltaOn  = IMM.pVS_DeltaOn;
            newIMM.pVS_Toff     = IMM.pVS_Toff;
            newIMM.pEnsureOneActive = IMM.pEnsureOneActive;

            newIMM.pFastMode       = IMM.pFastMode;
            newIMM.pUseEntropyCtrl = IMM.pUseEntropyCtrl;
            newIMM.pEntropyTarget  = IMM.pEntropyTarget;
            newIMM.pSkipInactive   = IMM.pSkipInactive;
            newIMM.pInactiveTicks  = IMM.pInactiveTicks;
            newIMM.pModelKind      = IMM.pModelKind;

            % New params
            newIMM.pRobustKind      = IMM.pRobustKind;
            newIMM.pNu              = IMM.pNu;
            newIMM.pNu_Jam          = IMM.pNu_Jam;
            newIMM.pHuberK          = IMM.pHuberK;
            newIMM.pHuberK_Jam      = IMM.pHuberK_Jam;
            newIMM.pJAM_CUSUM       = IMM.pJAM_CUSUM;
            newIMM.pJAM_ThreshA     = IMM.pJAM_ThreshA;
            newIMM.pJAM_Drift       = IMM.pJAM_Drift;
            newIMM.pJamActive       = IMM.pJamActive;

            newIMM.pAdaptAlphaMax   = IMM.pAdaptAlphaMax;
            newIMM.pAdaptGLRGain    = IMM.pAdaptGLRGain;
            newIMM.pAdaptEntGain    = IMM.pAdaptEntGain;
            newIMM.pAdaptWinnerBias = IMM.pAdaptWinnerBias;
            newIMM.pBoostToCA       = IMM.pBoostToCA;
            newIMM.pBoostToCV       = IMM.pBoostToCV;

            newIMM.pEpsSafeWTA      = IMM.pEpsSafeWTA;
            newIMM.pEpsInactiveMass = IMM.pEpsInactiveMass;
            newIMM.pTopKMax         = IMM.pTopKMax;

            newIMM.pTPM_Base        = IMM.pTPM_Base;
            newIMM.pTPM_Current     = IMM.pTPM_Current;
        end

        %------------------------------ Retrodiction stubs -----------------
        function [retroState, retroCov, success] = retrodict(obj, dt)
            coder.internal.errorIf(numel(obj) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'retrodict');
            validateRetrodict(obj, dt);
            [~, tk, tkappa, success, ~, ~] = getRetrodictQuantities(obj, dt);
            if ~success
                retroState = obj.State;
                retroCov   = obj.StateCovariance;
                obj.pLastRetrodictionDT(1) = 0;
                return
            end
            for i = coder.unroll(1:obj.pNumModels)
                retrodict(obj.pTrackingFilters{i}, dt);
            end
            retroState = obj.State;
            retroCov   = obj.StateCovariance;
            obj.pLastRetrodictionDT(1) = tk - tkappa;
        end

        function [x_retroCorr, P_retroCorr] = retroCorrect(obj, z, varargin)
            coder.internal.errorIf(numel(obj) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'retroCorrect');
            coder.internal.assert(obj.pWasRetrodicted, ...
                'shared_tracking:OOSMFilter:MustCallThisMethodBeforeThatMethod', ...
                'retrodict', 'retroCorrect');
            [~, distk] = fetchDistributionByTime(obj, max(obj.pCorrectionTimestamps));
            for kMdl = coder.unroll(1:obj.pNumModels)
                retroCorrect(obj.pTrackingFilters{kMdl}, z, distk.FilterDistributions{kMdl});
            end
            x_retroCorr = obj.State;
            P_retroCorr = obj.StateCovariance;
        end

        function [x_retroCorr, P_retroCorr] = retroCorrectJPDA(obj, z, jpda, varargin)
            coder.internal.assert(obj.pWasRetrodicted, ...
                'shared_tracking:OOSMFilter:MustCallThisMethodBeforeThatMethod', ...
                'retrodict', 'retroCorrectJPDA');
            [~, distk] = fetchDistributionByTime(obj, max(obj.pCorrectionTimestamps));
            for kMdl = coder.unroll(1:obj.pNumModels)
                retroCorrectJPDA(obj.pTrackingFilters{kMdl}, z, jpda, distk.FilterDistributions{kMdl});
            end
            x_retroCorr = obj.State;
            P_retroCorr = obj.StateCovariance;
        end

        function tunableProps = tunableProperties(obj)
            s1 = transProbTunableProps(obj);
            s2 = modelProbTunableProps(obj);
            f  = cell(obj.pNumModels,1);
            for i = 1:obj.pNumModels
                f{i} = tunableProperties(obj.pTrackingFilters{i});
            end
            tunableProps = tunableFilterProperties("trackingSAFEIMM", [s1; s2], f);
        end

        function setTunedProperties(obj, s)
            validateattributes(s, {'struct'}, {'nonempty'}, 'setTunedProperties', 'S');
            obj.TransitionProbabilities = s.TransitionProbabilities ./ sum(s.TransitionProbabilities, 2);
            obj.ModelProbabilities      = s.ModelProbabilities ./ sum(s.ModelProbabilities);
            for i = coder.unroll(1:obj.pNumModels)
                setTunedProperties(obj.pTrackingFilters{i}, s.TrackingFilters{i});
            end
        end
    end

    %==================== Access for base/contains/test ====================
    methods (Access = ...
            {?matlabshared.tracking.internal.RetrodictionFilter, ...
             ?matlabshared.tracking.internal.AbstractContainsFilters, ...
             ?matlab.unittest.TestCase})

        function retrodictionImpl(obj, dt, tk, tkappa, ~, ~, ~)
            for i = coder.unroll(1:obj.pNumModels)
                retrodict(obj.pTrackingFilters{i}, dt);
            end
            obj.pWasRetrodicted = true;
            obj.pLastRetrodictionDT(1) = tk - tkappa;
        end

        function retroCorrectImpl(obj, z, distk, varargin)
            likelihoods = zeros(obj.pNumModels,1,'like',obj.pModelProbabilities);
            coder.unroll();
            for kMdl = 1:obj.pNumModels
                fk  = obj.pTrackingFilters{kMdl};
                fdk = distk.FilterDistributions{kMdl};
                % Keep retro pipeline Gaussian; robust path is in online correct()
                likelihoods(kMdl) = likelihood(fk, z);
                retroCorrectImpl(fk, z, fdk);
            end
            cbar      = distk.ModelProbabilities(:);
            % Fast priors: Π_k' * cbar (use current Π)
            mdl_probs = obj.pTPM_Current.' * cbar;
            % Apply likelihoods
            mdl_probs = mdl_probs .* likelihoods;
            mdl_probs = mdl_probs / sum(mdl_probs);

            obj.pModelProbabilities = mdl_probs;
            [obj.pState, obj.pStateCovariance] = obj.weightedCombine(obj.pModelProbabilities);

            obj.pWasRetrodicted        = false;
            obj.pLastRetrodictionDT(1) = 0;
        end

        function [x_corr, P_corr] = retroCorrectJPDAImpl(obj, z, jpda, distk, varargin)
            M = size(z,2);
            likelihoods = zeros(obj.pNumModels, M, 'like', obj.pModelProbabilities);
            for kMdl = coder.unroll(1:obj.pNumModels)
                fk  = obj.pTrackingFilters{kMdl};
                fdk = distk.FilterDistributions{kMdl};
                for m = 1:M
                    likelihoods(kMdl,m) = likelihood(fk, z(:,m));
                end
                retroCorrectJPDAImpl(fk, z, jpda, fdk);
            end
            updateModelProbJPDA(obj, likelihoods, jpda, M);
            [x_corr, P_corr] = obj.weightedCombine(obj.pModelProbabilities);
            obj.pState = x_corr; obj.pStateCovariance = P_corr;

            obj.pWasRetrodicted        = false;
            obj.pLastRetrodictionDT(1) = 0;
        end
    end

    methods (Access = ...
            {?matlabshared.tracking.internal.AbstractTrackingFilter, ...
             ?matlabshared.tracking.internal.AbstractContainsFilters, ...
             ?matlab.unittest.TestCase})
        function sync(IMM, IMM2)
            classToExpect = 'trackingSAFEIMM';
            validateattributes(IMM2, {classToExpect}, {'scalar'}, 'trackingSAFEIMM');
            for i = coder.unroll(1:IMM.pNumModels)
                sync(IMM.pTrackingFilters{i}, IMM2.pTrackingFilters{i});
            end
            syncRetroFilter(IMM, IMM2);

            IMM.ModelConversionFcn     = IMM2.ModelConversionFcn;
            IMM.pModelProbabilities    = IMM2.pModelProbabilities;
            IMM.TransitionProbabilities = IMM2.TransitionProbabilities;
            IMM.MeasurementNoise       = IMM2.MeasurementNoise;

            [IMM.pState, IMM.pStateCovariance] = IMM.weightedCombine(IMM.pModelProbabilities);

            % sync new matrices
            IMM.pTPM_Base    = IMM2.pTPM_Base;
            IMM.pTPM_Current = IMM2.pTPM_Current;
        end

        function nullify(IMM)
            coder.internal.errorIf(numel(IMM) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', 'nullify');
            for kMdl = coder.unroll(1:IMM.pNumModels)
                nullify(IMM.pTrackingFilters{kMdl});
            end
            createHistory(IMM);
            [IMM.pState, IMM.pStateCovariance] = IMM.weightedCombine(IMM.pModelProbabilities);
        end

        function names = modelName(IMM), names = IMM.ModelNames; end

        function [stm, mm] = models(IMM, dt)
            [stm, mm] = models(IMM.pTrackingFilters{1}, dt);
        end

        function tf = supportsVarsizeMeasurements(obj)
            tf = true;
            for i = coder.unroll(1:numel(obj.pNumModels))
                tf = tf && supportsVarsizeMeasurements(obj.pTrackingFilters{i});
                if ~tf, return; end
            end
        end
    end

    %=========================== Protected utils ===========================
    methods (Access = protected)
        function filter = defaultTrackingFilters(~)
            % Minimal placeholders; your initFcn should override.
            filter = { ...
                trackingKF, ...
                trackingKF, ...
                trackingEKF };
        end

        function value = defaultModelConversionFcn(~)
            value = @switchimm;
        end

        function validateCorrectJPDAInputs(filter, z, jpda, fcname, varargin)
            coder.internal.errorIf(numel(filter) > 1, 'fusion:trackingSAFEIMM:NonScalarFilter', 'trackingSAFEIMM', fcname);
            validateattributes(z, {'double','single'}, {'real','finite','nonsparse','2d'}, fcname, 'z');
            ONE      = matlabshared.tracking.internal.fusion.codegen.StrictSingleUtilities.IntOne();
            numMeas  = matlabshared.tracking.internal.fusion.codegen.StrictSingleUtilities.IntIndex(size(z,2));
            expectedNumel = numMeas + ONE;
            validateattributes(jpda, {'double','single'}, ...
                {'real','finite','nonsparse','nonnegative','vector','numel',expectedNumel,'<=',1}, fcname, 'jpda');
            classToUse = class(filter.pTrackingFilters{1}.State);
            coder.internal.errorIf(abs(sum(jpda) - 1) > sqrt(eps(classToUse)), ...
                'shared_tracking:ExtendedKalmanFilter:InvalidJPDAInput', 'jpda')
        end

        function mdl_probs = updateModelProbJPDA(IMM, likelihoods, jpda, M)
            cbar     = IMM.pModelProbabilities;
            jpda_col = jpda(:);
            beta     = jpda_col(1:M);
            beta_not = jpda_col(M+1);
            jpdaLikelihoods = bsxfun(@plus, likelihoods*beta, beta_not);
            mdl_probs       = cbar(:) .* jpdaLikelihoods;
            mdl_probs       = mdl_probs/sum(mdl_probs);
            IMM.pModelProbabilities = mdl_probs;
        end
    end

    %===================== OOSM history support ============================
    methods (Access = {?matlabshared.tracking.internal.OOSMFilter, ?matlab.unittest.TestCase})
        function createHistory(obj)
            if obj.MaxNumOOSMSteps > 0
                dist      = getDistribution(obj);
                classToUse = class(dist.FilterDistributions{1}.State);
                obj.pPredictionDeltaTFromLastCorrection = zeros(1,1,classToUse);
                obj.pCorrectionTimestamps               = zeros(1, obj.MaxNumOOSMSteps, classToUse);
                obj.pCorrectionDistributions            = repmat(dist, 1, obj.MaxNumOOSMSteps);
                obj.pLastRetrodictionDT                 = zeros(1,1,classToUse);
            end
        end
    end

    %========================== Private helpers ===========================
    methods (Access = private)

        function oh = oneHot(IMM, idx)
            % Return e_i as a column vector in R^{pNumModels}
            % Uses the same numeric class as pModelProbabilities (codegen safe)
            L  = IMM.pNumModels;
            oh = zeros(L,1,'like',IMM.pModelProbabilities);
            idx = max(1, min(double(idx), double(L))); % clamp for safety
            oh(idx) = cast(1,'like',oh);
        end

        function updateActiveMask(IMM, mdl_probs, widx)
            % Update per-model active mask using hysteresis counters.
            % Turn OFF if prob stays below pVS_DeltaOff for pVS_Toff steps.
            % Turn ON immediately if prob rises above pVS_DeltaOn.
            % Always keep winner active; ensure at least one model active.

            % increment low-probability counters / deactivate if needed
            for i = 1:IMM.pNumModels
                if mdl_probs(i) < IMM.pVS_DeltaOff
                    IMM.pLowCount(i) = IMM.pLowCount(i) + 1;
                else
                    IMM.pLowCount(i) = uint16(0);
                end
                if IMM.pLowCount(i) >= IMM.pVS_Toff
                    IMM.pActiveMask(i) = false;
                end
            end

            % keep the current winner active
            IMM.pActiveMask(widx) = true;

            % reactivate any model whose prob exceeds the ON threshold
            for i = 1:IMM.pNumModels
                if mdl_probs(i) > IMM.pVS_DeltaOn
                    IMM.pActiveMask(i) = true;
                end
            end

            % safety: ensure at least one model is active
            if IMM.pEnsureOneActive && ~any(IMM.pActiveMask)
                IMM.pActiveMask(widx) = true;
            end
        end
        % Priors via TPM: c <- Π_k' * c
        function advancePriors(IMM)
            IMM.pModelProbabilities(:) = IMM.pTPM_Current.' * IMM.pModelProbabilities(:);
        end

        % After predict, provide a neutral readout (use last winner)
        function [x_pred, P_pred] = readoutAfterPredict(IMM)
            w = double(IMM.pBC_LastWinner); w = max(1, min(w, IMM.pNumModels));
            [x_pred, P_pred] = IMM.readFromModel(w);
            IMM.pState = x_pred; IMM.pStateCovariance = P_pred;
        end

        % Convert (xi,Pi) from model i -> model j space using ModelConversionFcn only if needed
        function [xr, Pr] = convertStateCov(IMM, fromIdx, xi, Pi, toIdx)
            refF = IMM.pTrackingFilters{toIdx};
            if numel(xi) == numel(refF.State)
                xr = xi; Pr = Pi; return
            end
            Xzeros = zeros(size(refF.State),           'like', refF.State);
            Pzeros = zeros(size(refF.StateCovariance), 'like', refF.StateCovariance);
            xr = IMM.ModelConversionFcn(IMM.ModelNames{fromIdx}, xi, IMM.ModelNames{toIdx}, Xzeros);
            if nargin(IMM.ModelConversionFcn) == 5
                Pr = IMM.ModelConversionFcn(IMM.ModelNames{fromIdx}, Pi, IMM.ModelNames{toIdx}, Pzeros, xi);
            else
                Pr = IMM.ModelConversionFcn(IMM.ModelNames{fromIdx}, Pi, IMM.ModelNames{toIdx}, Pzeros);
            end
        end

        % Weighted combination over a subset (moment matching)
        function [Xout, Pout] = weightedCombine(IMM, weights, keepIdx)
            if nargin < 3, keepIdx = 1:IMM.pNumModels; end

            w = weights(:);
            mask = false(IMM.pNumModels,1); mask(keepIdx) = true;
            for i=1:numel(w), if ~mask(i), w(i)=0; end, end
            s = sum(w);
            if s <= 0, w = ones(IMM.pNumModels,1,'like',w)/IMM.pNumModels; else, w = w/s; end

            refIdx = 1;
            refF   = IMM.pTrackingFilters{refIdx};

            Xout = zeros(size(refF.State),           'like', refF.State);
            Pout = zeros(size(refF.StateCovariance), 'like', refF.StateCovariance);

            % First moment
            for i = 1:IMM.pNumModels
                wi = w(i); if wi == 0, continue; end
                fi = IMM.pTrackingFilters{i};
                [xi_ref, ~] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, refIdx);
                Xout = Xout + wi * xi_ref;
            end

            % Second moment
            for i = 1:IMM.pNumModels
                wi = w(i); if wi == 0, continue; end
                fi = IMM.pTrackingFilters{i};
                [xi_ref, Pi_ref] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, refIdx);
                dx   = xi_ref - Xout;
                Pout = Pout + wi * (Pi_ref + dx*dx.');
            end
        end

        function [x, P] = readFromModel(IMM, idx)
            x = IMM.pTrackingFilters{idx}.State;
            P = IMM.pTrackingFilters{idx}.StateCovariance;
        end

        % Cheap posterior sharpness: pmax>=0.75 and pmax/p2 >= 3
        function tf = posteriorIsSharp(~, pmax, p2)
            if p2==0, r = inf; else, r = pmax/p2; end
            tf = (pmax >= 0.75) && (r >= 3);
        end

        % Return best and second best probabilities and values
        function [pmax, imax, p2] = bestTwoWithValues(~, p)
            imax = 1; i2 = 1; a = -inf; b = -inf;
            for k = 1:numel(p)
                v = p(k);
                if v > a
                    b = a; i2 = imax;
                    a = v; imax = k;
                elseif v > b
                    b = v; i2 = k;
                end
            end
            pmax = a; p2 = b;
        end

        % Minimum Mahalanobis^2 separation between winner and rivals (Cholesky)
        function m2 = minMah2ToRivals(IMM, widx)
            xw = IMM.pTrackingFilters{widx}.State;
            Pw = IMM.pTrackingFilters{widx}.StateCovariance;

            m2 = inf;
            for i = 1:IMM.pNumModels
                if i == widx || ~IMM.pActiveMask(i), continue; end
                fi = IMM.pTrackingFilters{i};
                [xiw, Piw] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, widx);
                Pbar = 0.5*(Pw + Piw);
                % Robust symmetric PSD & jitter
                Pbar = (Pbar+Pbar.')*0.5;
                ridge = single(1e-6);
                Pbar  = Pbar + ridge*eye(size(Pbar), 'like', Pbar);
                % Cholesky
                R = chol(Pbar,'lower');
                dx = xiw - xw;
                y  = R \ dx;
                z  = R.' \ y;
                m2_i = real(dx' * z);
                if m2_i < m2, m2 = m2_i; end
            end
            if isinf(m2), m2 = 0; end
        end

        function B = epsSafeBoundWTA(IMM, widx, weights)
        % ε-safe WTA bound (Option A: all rivals included)
        %   B = t * sqrt( trace(Pbar) * avg_d2 )
        % where
        %   t       = 1 - w_winner
        %   Pbar    = 0.5 * ( Pw + sum_{i≠w} (w_i/t) * P_i→w )
        %   avg_d2  = (1/t) * sum_{i≠w} w_i * || μ_i→w - μ_w ||^2_{Pbar^{-1}}
        %
        % All means/covs are mapped into the winner space (via convertStateCov).
        % Codegen-safe: no cells/try-catch; fixed-size ops; robust Cholesky with ridge.

            % ----- normalize weights -----
            M = IMM.pNumModels;
            w = weights(:);
            s = sum(w);
            if s <= 0
                w = ones(M,1,'like',w) / M;
            else
                w = w / s;
            end

            % ----- clamp winner index -----
            widx = max(1, min(double(widx), double(M)));

            % ----- winner state/cov -----
            fw  = IMM.pTrackingFilters{widx};
            muw = fw.State;
            Pw  = fw.StateCovariance;

            % ----- tail mass -----
            tw = w(widx);
            t  = max(0, 1 - tw);
            if t <= eps('like', t)
                B = cast(0, 'like', t);
                return
            end

            % ----- first pass: build Σ (w_i/t) * P_i→w over ALL rivals -----
            sumPi = zeros(size(Pw), 'like', Pw);
            for i = 1:M
                if i == widx, continue; end
                fi = IMM.pTrackingFilters{i};
                [~, Pi_w] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, widx);
                sumPi = sumPi + (w(i)/t) * Pi_w;
            end

            % ----- Pbar (symmetrized) -----
            Pbar = 0.5 * (Pw + sumPi);
            Pbar = 0.5 * (Pbar + Pbar.');  % numerically symmetrize
            ridge = cast(1e-6, 'like', Pbar);
            Ibar  = eye(size(Pbar), 'like', Pbar);
            trP   = trace(Pbar);
            if trP <= eps('like',trP), trP = eps('like',trP); end

            % ----- rival-set hash (guards cache reuse when set changes) -----
            % (All rivals included: hash depends only on indices 1..M except winner.)
            h = uint32(0);
            for i = 1:M
                if i ~= widx
                    h = h * uint32(1664525) + uint32(i);
                end
            end

            % ----- decide reuse vs recompute BEFORE using R -----
            reuseOK = IMM.pEPS_CacheValid && ...
                      (IMM.pEPS_LastWinner == uint32(widx)) && ...
                      (IMM.pEPS_LastHash   == h);

            tolRel = cast(0.25,'like',trP);   % 25% relative change allowed for reuse
            if IMM.pEPS_CacheValid && (IMM.pEPS_LastWinner == uint32(widx))
                rel = abs(trP - IMM.pEPS_trPbar) / max(IMM.pEPS_trPbar, eps('like',trP));
            else
                rel = inf;
            end

            % debug and audit
            ridgemult = uint16(1);
            [R,p] = chol(Pbar + ridge*Ibar, 'lower');
            if p > 0
                ridgemult = uint16(10);
                [R,p] = chol(Pbar + 10*ridge*Ibar, 'lower');
                if p > 0
                    ridgemult = uint16(100);
                    [R,p] = chol(Pbar + 100*ridge*Ibar, 'lower');
                end
            end

            % cache (your existing cache writes) ...
            IMM.pEPS_R_chol     = R;
            IMM.pEPS_trPbar     = cast(trP, 'like', IMM.pEPS_trPbar);
            IMM.pEPS_LastWinner = uint32(widx);
            IMM.pEPS_CacheValid = true;

            % debug (only if enabled; otherwise free)
            if IMM.pDBG_Enable
                IMM.pDBG_RidgeMult = ridgemult;
            end



            if reuseOK && (rel <= tolRel)
                R = IMM.pEPS_R_chol;
            else
                % robust Cholesky with escalating ridge (codegen-safe)
                [R,p] = chol(Pbar + ridge*Ibar, 'lower');
                if p > 0
                    [R,p] = chol(Pbar + 10*ridge*Ibar, 'lower');
                    if p > 0
                        [R,p] = chol(Pbar + 100*ridge*Ibar, 'lower');
                    end
                end
                % update cache
                IMM.pEPS_R_chol     = R;
                IMM.pEPS_trPbar     = cast(trP, 'like', IMM.pEPS_trPbar);
                IMM.pEPS_LastWinner = uint32(widx);
                IMM.pEPS_LastHash   = h;
                IMM.pEPS_CacheValid = true;
            end

            % ----- second pass: average Mahalanobis^2 under Pbar over ALL rivals -----
            acc_d2 = zeros(1,1,'like',t);
            for i = 1:M
                if i == widx, continue; end
                fi = IMM.pTrackingFilters{i};
                [mui_w, ~] = IMM.convertStateCov(i, fi.State, fi.StateCovariance, widx);
                dx = mui_w - muw;
                y  = R \ dx;
                z  = R.' \ y;
                d2 = max(real(dx' * z), 0);
                acc_d2 = acc_d2 + w(i) * d2;
            end
            avg_d2 = acc_d2 / t;

            % After computing trP and avg_d2, just before final B:
            IMM.pDBG_trPbar = single(trP);
            IMM.pDBG_avgd2  = single(avg_d2);


            % ----- final bound -----
            B = cast(t, 'like', trP) * sqrt(trP * avg_d2);
            if ~isfinite(B), B = realmax('like', trP); end
        end

        %=================== NEW: ε-safe active mask =======================
        function updateActiveMaskSafe(IMM, mdl_probs, widx)
            % First do the usual hysteresis-based on/off
            for i = 1:IMM.pNumModels
                if mdl_probs(i) < IMM.pVS_DeltaOff
                    IMM.pLowCount(i) = IMM.pLowCount(i) + 1;
                else
                    IMM.pLowCount(i) = uint16(0);
                end
                if IMM.pLowCount(i) >= IMM.pVS_Toff
                    IMM.pActiveMask(i) = false;
                end
            end
            IMM.pActiveMask(widx) = true; % keep winner
            for i = 1:IMM.pNumModels
                if mdl_probs(i) > IMM.pVS_DeltaOn
                    IMM.pActiveMask(i) = true;
                end
            end
            if IMM.pEnsureOneActive && ~any(IMM.pActiveMask)
                IMM.pActiveMask(widx) = true;
            end

            % ε-safe mass cap: ensure inactive mass ≤ pEpsInactiveMass
            inactMass = sum(mdl_probs(~IMM.pActiveMask));
            if inactMass > IMM.pEpsInactiveMass
                % Reactivate highest-prob inactives until mass bound satisfied
                idxInact = find(~IMM.pActiveMask);
                % simple insertion sort by prob
                for ii = 1:numel(idxInact)-1
                    for jj = ii+1:numel(idxInact)
                        if mdl_probs(idxInact(jj)) > mdl_probs(idxInact(ii))
                            tmp = idxInact(ii); idxInact(ii) = idxInact(jj); idxInact(jj) = tmp;
                        end
                    end
                end
                k = 1;
                while inactMass > IMM.pEpsInactiveMass && k <= numel(idxInact)
                    IMM.pActiveMask(idxInact(k)) = true;
                    inactMass = sum(mdl_probs(~IMM.pActiveMask));
                    k = k + 1;
                end
            end
        end

        %=================== NEW: robust likelihoods =======================
        function lh = robustLikelihoodFromNIS(IMM, d2, mDim)
            % mDim = measurement dimension
            nu     = IMM.pNu;
            huberK = IMM.pHuberK;

            if IMM.pJamActive
                % Harden under jam
                if IMM.pRobustKind == 1
                    nu = IMM.pNu_Jam;
                elseif IMM.pRobustKind == 2
                    huberK = IMM.pHuberK_Jam;
                end
            end

            switch IMM.pRobustKind
                case uint8(1) % Student-t (unnormalized, sufficient for ratios)
                    % lhood ∝ (1 + d2/nu)^(-(nu + m)/2)
                    lh = (1 + d2/max(nu,eps('single')))^(-0.5*(nu + mDim));
                case uint8(2) % Huber pseudo-likelihood: exp(-rho/2)
                    r = sqrt(max(d2,0));
                    k = huberK;
                    if r <= k
                        rho = r*r;
                    else
                        rho = 2*k*r - k*k;
                    end
                    lh = exp(-0.5*rho);
                otherwise     % Gaussian fallback using d2 only
                    % lhood ∝ exp(-0.5*d2)
                    lh = exp(-0.5*max(d2,0));
            end
            if ~isfinite(lh), lh = realmin('single'); end
            lh = max(lh, realmin('single'));
        end

        function d2 = safeDistance(IMM, fk, z, varargin)
        % Prefer distance(fk,z,measParams) if available; fall back to likelihood.

            if coder.target('MATLAB')
                % Try with meas params first (some filters require them)
                try
                    d2 = distance(fk, z, varargin{:});
                catch
                    try
                        d2 = distance(fk, z);
                    catch
                        % Fallback: map likelihood to pseudo-distance
                        try
                            lg = likelihood(fk, z, varargin{:});
                            lg = max(lg, realmin('like', single(1)));
                            d2 = -2*log(lg);
                        catch
                            d2 = single(0);
                        end
                    end
                end
            else
                % Codegen path: assume supported tracking filters implement distance()
                d2 = distance(fk, z);
            end
            if ~isfinite(d2), d2 = single(0); end
        end



        %=================== NEW: Adaptive Π_k =============================
        function updateAdaptiveTPM(IMM)
            % Compute alpha from GLR and entropy of current model probs
            mu = IMM.pModelProbabilities(:);
            mu = mu / sum(mu);
            % Shannon entropy (nats)
            H  = -sum( max(mu,realmin('single')) .* log(max(mu,realmin('single'))) );
            % Normalize entropy to [0,1] by dividing by log(M)
            Hn = H / log(max(single(IMM.pNumModels),1));

            % GLR-derived component (larger GLR_CUSUM -> higher alpha)
            a_glr = min(IMM.pAdaptGLRGain * max(IMM.pGLR_CUSUM,0), 1);
            % Entropy-derived (low entropy -> less mixing; high entropy -> more)
            a_ent = min(IMM.pAdaptEntGain * max(Hn,0), 1);

            alpha = min(IMM.pAdaptAlphaMax, a_glr + a_ent);

            % Build a boost matrix that favors CA-like under maneuvers
            B = IMM.pTPM_Base; % start from base topology
            % Small bias to move out of current winner if evidence is strong
            if any(IMM.pActiveMask)
                w = double(IMM.pBC_LastWinner);
                for j = 1:IMM.pNumModels
                    if j == w
                        % reduce self-stay a bit under high alpha
                        B(j,j) = max(B(j,j) - IMM.pAdaptWinnerBias*alpha, 0);
                    end
                end
            end
            % Row-wise boosts: CV rows push more prob to CA columns
            for r = 1:IMM.pNumModels
                if IMM.pModelKind(r) == uint8(1) % CV row
                    add = IMM.pBoostToCA * alpha;
                    % distribute 'add' equally to CA columns
                    idxCA = (IMM.pModelKind == uint8(2));
                    if any(idxCA)
                        share = add / sum(idxCA);
                        B(r,idxCA) = B(r,idxCA) + share;
                        % subtract from self to keep row sum ~ 1 before renorm
                        B(r,r) = max(B(r,r) - add, 0);
                    end
                else % CA row: when calm (low alpha), nudge towards CV
                    add = IMM.pBoostToCV * max(0, IMM.pAdaptAlphaMax - alpha);
                    idxCV = (IMM.pModelKind == uint8(1));
                    if any(idxCV)
                        share = add / sum(idxCV);
                        B(r,idxCV) = B(r,idxCV) + share;
                        B(r,r) = max(B(r,r) - add, 0);
                    end
                end
            end

            % Blend and renormalize rows strictly
            T = (1 - alpha)*IMM.pTPM_Base + alpha*B;
            for r = 1:IMM.pNumModels
                sr = sum(T(r,:));
                if sr <= 0
                    T(r,:) = 1/IMM.pNumModels;
                else
                    T(r,:) = T(r,:) / sr;
                end
            end
            IMM.pTPM_Current = T;
        end

        %=================== NEW: top-K selector ===========================
        function keepIdx = selectTopK(~, p, K)
            % Return indices of top-K probabilities (K>=1)
            K = max(1, min(double(K), numel(p)));
            keepIdx = zeros(1,K,'uint32');
            used = false(size(p));
            for k = 1:K
                imax = 1; a = -inf;
                for i = 1:numel(p)
                    if ~used(i)
                        v = p(i);
                        if v > a
                            a = v; imax = i;
                        end
                    end
                end
                keepIdx(k) = uint32(imax);
                used(imax) = true;
            end
        end

        % ================= BESTCHOICE-IMM helpers (kept) ==================
        function [vel, yaw, hasVel] = bcExtractVelYaw(IMM, state)
            vel    = zeros(3,1,'like',state);
            yaw    = single(0);
            hasVel = false;
            if numel(state) >= 6
                vx  = state(2);
                vy  = state(4);
                vz  = state(6);
                vel = [vx; vy; vz];
                spd = hypot(vx, vy);
                if spd > 1e-3
                    yaw    = atan2(vy, vx);
                    hasVel = true;
                end
            end
        end

        function F = bcMakeFcv(~, dt, likeX)
            dt = cast(dt, 'like', likeX);
            F  = blkdiag([1 dt;0 1],[1 dt;0 1],[1 dt;0 1]);
        end

        function Q = bcMakeQcv(IMM, dt, likeX)
            dt = cast(dt, 'like', likeX);
            q  = cast(IMM.pCVq, 'like', likeX);
            Qa = q * [dt^4/4, dt^3/2; dt^3/2, dt^2];
            Q  = blkdiag(Qa, Qa, Qa);
        end

        function F = bcMakeFca(~, dt, likeX)
            dt = cast(dt, 'like', likeX);
            F3 = [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];
            F  = blkdiag(F3, F3, F3);
        end

        function Q = bcMakeQca(IMM, dt, likeX)
            dt = cast(dt, 'like', likeX);
            q  = cast(IMM.pCAq, 'like', likeX);
            Qa = q * [dt^5/20, dt^4/8, dt^3/6; dt^4/8, dt^3/3, dt^2/2; dt^3/6, dt^2/2, dt];
            Q  = blkdiag(Qa, Qa, Qa);
        end

        %====================== JAMMING CUSUM ==============================
        function updateJammingCUSUM(IMM, avgNIS, mDim)
            % Expected NIS ≈ mDim under Gaussian assumptions.
            mRef = max(single(1), cast(mDim, 'like', avgNIS));
            incr = max(single(0), avgNIS/max(mRef,1) - 1.0) - IMM.pJAM_Drift;

            IMM.pJAM_CUSUM = max(single(0), IMM.pJAM_CUSUM + incr);
            IMM.pJamActive = IMM.pJAM_CUSUM >= IMM.pJAM_ThreshA;

            if ~IMM.pJamActive
                % small decay toward zero
                IMM.pJAM_CUSUM = max(single(0), IMM.pJAM_CUSUM - 0.1*IMM.pJAM_Drift);
            end
        end

    end

    %=========================== Custom display ============================
    methods (Access = 'protected')
        function propGroups = getPropertyGroups(obj)
            propGroups = [ ...
                matlab.mixin.util.PropertyGroup({'State','StateCovariance'}, ''), ...
                matlab.mixin.util.PropertyGroup({'TrackingFilters','ModelNames','HasMeasurementWrapping','MeasurementNoise'}, ''), ...
                matlab.mixin.util.PropertyGroup({'ModelConversionFcn','TransitionProbabilities','ModelProbabilities'}, '') ...
            ];
            oosmGroup    = getPropertyGroups@matlabshared.tracking.internal.RetrodictionFilter(obj);
            smoothGroups = getPropertyGroups@fusion.internal.IMMSmoother(obj);
            propGroups   = [propGroups; oosmGroup; smoothGroups];
        end
    end

    %=============================== Save/load =============================
    methods (Access = protected)
        function [sobj] = saveobj(IMM)
            sIMM = struct( ...
                'State',               IMM.State, ...
                'StateCovariance',     IMM.StateCovariance, ...
                'ModelConversionFcn',  IMM.ModelConversionFcn, ...
                'TransitionProbabilities', IMM.TransitionProbabilities, ...
                'MeasurementNoise',    IMM.MeasurementNoise, ...
                'ModelProbabilities',  IMM.ModelProbabilities, ...
                'pState',              IMM.pState, ...
                'pStateCovariance',    IMM.pStateCovariance, ...
                'pTransitionProbabilities', IMM.pTransitionProbabilities, ...
                'pTPM_Base',           IMM.pTPM_Base, ...
                'pTPM_Current',        IMM.pTPM_Current, ...
                'pMeasurementNoise',   IMM.pMeasurementNoise, ...
                'pModelProbabilities', IMM.pModelProbabilities, ...
                'pNumModels',          IMM.pNumModels, ...
                'TrackingFilters',     {IMM.pTrackingFilters}, ...
                'pLastRetrodictionDT', IMM.pLastRetrodictionDT, ...
                'ModelNames',          {IMM.ModelNames}, ...
                'pActiveMask',         IMM.pActiveMask, ...
                'pLowCount',           IMM.pLowCount, ...
                'pGLR_CUSUM',          IMM.pGLR_CUSUM, ...
                'pJAM_CUSUM',          IMM.pJAM_CUSUM, ...
                'pJamActive',          IMM.pJamActive, ...
                'pFastMode',           IMM.pFastMode, ...
                'pUseEntropyCtrl',     IMM.pUseEntropyCtrl, ...
                'pEntropyTarget',      IMM.pEntropyTarget, ...
                'pSkipInactive',       IMM.pSkipInactive, ...
                'pInactiveTicks',      IMM.pInactiveTicks, ...
                'pModelKind',          IMM.pModelKind, ...
                'pRobustKind',         IMM.pRobustKind, ...
                'pNu',                 IMM.pNu, ...
                'pNu_Jam',             IMM.pNu_Jam, ...
                'pHuberK',             IMM.pHuberK, ...
                'pHuberK_Jam',         IMM.pHuberK_Jam, ...
                'pAdaptAlphaMax',      IMM.pAdaptAlphaMax, ...
                'pAdaptGLRGain',       IMM.pAdaptGLRGain, ...
                'pAdaptEntGain',       IMM.pAdaptEntGain, ...
                'pAdaptWinnerBias',    IMM.pAdaptWinnerBias, ...
                'pBoostToCA',          IMM.pBoostToCA, ...
                'pBoostToCV',          IMM.pBoostToCV, ...
                'pEpsSafeWTA',         IMM.pEpsSafeWTA, ...
                'pEpsInactiveMass',    IMM.pEpsInactiveMass, ...
                'pTopKMax',            IMM.pTopKMax );
            sSmooth = saveSmootherProperties(IMM, sIMM);
            sobj    = saveRetroProperties(IMM, sSmooth);
        end

        function loadPrivateProtectedProperties(IMM, sobj)
            IMM.pState                  = sobj.pState;
            IMM.pStateCovariance        = sobj.pStateCovariance;
            IMM.ModelConversionFcn      = sobj.ModelConversionFcn;
            IMM.pTransitionProbabilities = sobj.pTransitionProbabilities;
            IMM.pTPM_Base               = sobj.pTPM_Base;
            IMM.pTPM_Current            = sobj.pTPM_Current;
            IMM.pMeasurementNoise       = sobj.pMeasurementNoise;
            IMM.pModelProbabilities     = sobj.pModelProbabilities;
            IMM.pNumModels              = sobj.pNumModels;
            if isfield(sobj, 'pLastRetrodictionDT')
                IMM.pLastRetrodictionDT = sobj.pLastRetrodictionDT;
            else
                IMM.pLastRetrodictionDT = zeros(1,1,'like',sobj.pState);
            end
            if isfield(sobj, 'ModelNames')
                IMM.ModelNames = sobj.ModelNames;
            else
                IMM.ModelNames = defaultModelNames(IMM.TrackingFilters);
            end
            if isfield(sobj, 'pActiveMask'), IMM.pActiveMask = sobj.pActiveMask; end
            if isfield(sobj, 'pLowCount'),   IMM.pLowCount   = sobj.pLowCount;   end
            if isfield(sobj, 'pGLR_CUSUM'),  IMM.pGLR_CUSUM  = sobj.pGLR_CUSUM;  end
            if isfield(sobj, 'pJAM_CUSUM'),  IMM.pJAM_CUSUM  = sobj.pJAM_CUSUM;  end
            if isfield(sobj, 'pFastMode'),       IMM.pFastMode       = sobj.pFastMode; end
            if isfield(sobj, 'pUseEntropyCtrl'), IMM.pUseEntropyCtrl = sobj.pUseEntropyCtrl; end
            if isfield(sobj, 'pEntropyTarget'),  IMM.pEntropyTarget  = sobj.pEntropyTarget;  end
            if isfield(sobj, 'pSkipInactive'),   IMM.pSkipInactive   = sobj.pSkipInactive;   end
            if isfield(sobj, 'pInactiveTicks'),  IMM.pInactiveTicks  = sobj.pInactiveTicks;  end
            if isfield(sobj, 'pModelKind'),      IMM.pModelKind      = sobj.pModelKind;      end

            % New
            if isfield(sobj,'pRobustKind'),      IMM.pRobustKind      = sobj.pRobustKind;      end
            if isfield(sobj,'pNu'),              IMM.pNu              = sobj.pNu;              end
            if isfield(sobj,'pNu_Jam'),          IMM.pNu_Jam          = sobj.pNu_Jam;          end
            if isfield(sobj,'pHuberK'),          IMM.pHuberK          = sobj.pHuberK;          end
            if isfield(sobj,'pHuberK_Jam'),      IMM.pHuberK_Jam      = sobj.pHuberK_Jam;      end
            if isfield(sobj,'pAdaptAlphaMax'),   IMM.pAdaptAlphaMax   = sobj.pAdaptAlphaMax;   end
            if isfield(sobj,'pAdaptGLRGain'),    IMM.pAdaptGLRGain    = sobj.pAdaptGLRGain;    end
            if isfield(sobj,'pAdaptEntGain'),    IMM.pAdaptEntGain    = sobj.pAdaptEntGain;    end
            if isfield(sobj,'pAdaptWinnerBias'), IMM.pAdaptWinnerBias = sobj.pAdaptWinnerBias; end
            if isfield(sobj,'pBoostToCA'),       IMM.pBoostToCA       = sobj.pBoostToCA;       end
            if isfield(sobj,'pBoostToCV'),       IMM.pBoostToCV       = sobj.pBoostToCV;       end
            if isfield(sobj,'pEpsSafeWTA'),      IMM.pEpsSafeWTA      = sobj.pEpsSafeWTA;      end
            if isfield(sobj,'pEpsInactiveMass'), IMM.pEpsInactiveMass = sobj.pEpsInactiveMass; end
            if isfield(sobj,'pTopKMax'),         IMM.pTopKMax         = sobj.pTopKMax;         end
        end

        function dist = getDistribution(obj)
            dist = distribution(obj);
        end

        function setMaxNumOOSMSteps(IMM, value)
            if value > 0
                isRetrodictionFilter = zeros(1, numel(IMM.pTrackingFilters), 'logical');
                coder.unroll();
                for i = 1:numel(IMM.pTrackingFilters)
                    isRetrodictionFilter(i) = isa(IMM.pTrackingFilters{i}, 'matlabshared.tracking.internal.RetrodictionFilter');
                end
                coder.internal.assert(all(isRetrodictionFilter), 'fusion:trackingSAFEIMM:expectedRetrodictionFilters');
            end
            setMaxNumOOSMSteps@matlabshared.tracking.internal.RetrodictionFilter(IMM, value);
        end
    end

    methods (Static = true)
        function retIMM = loadobj(sobj)
            if isscalar(sobj)
                retIMM = trackingSAFEIMM(sobj.TrackingFilters, sobj.ModelConversionFcn, sobj.TransitionProbabilities);
                loadPrivateProtectedProperties(retIMM, sobj);
            else
                retIMM = trackingSAFEIMM({sobj.TrackingFilters}, sobj(1).ModelConversionFcn, sobj(1).TransitionProbabilities);
                loadPrivateProtectedProperties(retIMM, sobj(1));
            end
            loadSmootherProperties(retIMM, sobj);
            loadRetroProperties(retIMM, sobj);
        end
    end

    methods
        function setMeasurementSizes(obj, measurementSize, measurementNoiseSize)
            validateattributes(measurementSize, {'numeric'}, {'real','positive','integer','scalar'}, 'trackingSAFEIMM');
            validateattributes(measurementNoiseSize, {'numeric'}, {'real','positive','integer','scalar'}, 'trackingSAFEIMM');
            for i = coder.unroll(1:obj.pNumModels)
                if ~isa(obj.pTrackingFilters{i}, 'trackingABF')
                    setMeasurementSizes(obj.pTrackingFilters{i}, measurementSize, measurementNoiseSize)
                end
            end
        end
    end

    methods (Static, Hidden)
        function props = matlabCodegenNontunableProperties(~)
            smootherProps = fusion.internal.IMMSmoother.matlabCodegenNontunableProperties;
            props         = {smootherProps{:}, 'pIsLinearKalmanFilter'};
        end

        function validateBounds(prop, r)
            switch prop
                case 'ModelProbabilities'
                    validatedBounds = all(r.LowerBound >= 0 & r.UpperBound <= 1);
                case 'TransitionProbabilities'
                    validatedBounds = all(r.LowerBound >= 0 & r.UpperBound <= 1);
                otherwise
                    return
            end
            assert(validatedBounds, ...
                message('fusion:tunableFilterProperties:ViolatesPropLimit', prop, mfilename));
        end
    end

    %==================== Parsing helpers (sim/codegen) ====================
    methods (Access = private)
        function [state, stateCov, trackingFilters, ...
                  transitionProbabilities, modelConversion, ...
                  modelProb, measNoise, modelNames] = parseInputs(IMM, varargin)

            firstNVIndex = matlabshared.tracking.internal.findFirstNVPair(varargin{:});
            coder.internal.errorIf(firstNVIndex > 4, 'fusion:trackingSAFEIMM:invalidInputsToConstructor', firstNVIndex - 1);

            switch firstNVIndex
                case 1
                    parsingParams = coder.internal.constantPreservingStruct( ...
                        'State', [], ...
                        'StateCovariance', [], ...
                        'TrackingFilters', [], ...
                        'ModelConversionFcn', IMM.defaultModelConversionFcn, ...
                        'TransitionProbabilities', [], ...
                        'ModelProbabilities', [], ...
                        'MeasurementNoise', []);
                case 2
                    parsingParams = coder.internal.constantPreservingStruct( ...
                        'State', [], ...
                        'StateCovariance', [], ...
                        'TrackingFilters', varargin{1}, ...
                        'ModelConversionFcn', IMM.defaultModelConversionFcn, ...
                        'TransitionProbabilities', [], ...
                        'ModelProbabilities', [], ...
                        'MeasurementNoise', []);
                case 3
                    parsingParams = coder.internal.constantPreservingStruct( ...
                        'State', [], ...
                        'StateCovariance', [], ...
                        'TrackingFilters', varargin{1}, ...
                        'ModelConversionFcn', varargin{2}, ...
                        'TransitionProbabilities', [], ...
                        'ModelProbabilities', [], ...
                        'MeasurementNoise', []);
                case 4
                    parsingParams = coder.internal.constantPreservingStruct( ...
                        'State', [], ...
                        'StateCovariance', [], ...
                        'TrackingFilters', varargin{1}, ...
                        'ModelConversionFcn', varargin{2}, ...
                        'TransitionProbabilities', varargin{3}, ...
                        'ModelProbabilities', [], ...
                        'MeasurementNoise', []);
            end

            if coder.target('MATLAB')
                [state, stateCov, trackingFiltersIn, transitionProbabilitiesIn, modelConversion, modelProb, measNoise, modelNames] ...
                    = parseInputsSimulation(parsingParams, varargin{firstNVIndex:end});
            else
                [state, stateCov, trackingFiltersIn, transitionProbabilitiesIn, modelConversion, modelProb, measNoise, modelNames] ...
                    = parseInputsCodegen(parsingParams, varargin{firstNVIndex:end});
            end

            if isempty(trackingFiltersIn)
                trackingFilters = IMM.defaultTrackingFilters;
            else
                trackingFilters = trackingFiltersIn;
            end

            if isempty(transitionProbabilitiesIn)
                transitionProbabilities = cast(IMM.defaultTransitionProbabilities, 'like', trackingFilters{1}.State);
            else
                transitionProbabilities = transitionProbabilitiesIn;
            end
        end
    end

    methods (Static)
        function names = defaultModelNames(trackingFilters)
            names = cell(1, numel(trackingFilters));
            coder.unroll();
            for i = 1:numel(trackingFilters)
                if ismethod(trackingFilters{i}, 'modelName')
                    names{i} = coder.const(trackingFilters{i}.modelName);
                else
                    names{i} = sprintf('model_%d', i);
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% Parse inputs for simulation
% -------------------------------------------------------------------------
function [state, stateCov, trackingFilters, transitionProbabilities, ...
          modelConversion, modelProb, measNoise, modelNames] = parseInputsSimulation(defaultParams, varargin)

parser = inputParser;
parser.addParameter('State',                defaultParams.State);
parser.addParameter('StateCovariance',      defaultParams.StateCovariance);
parser.addParameter('TrackingFilters',      defaultParams.TrackingFilters);
parser.addParameter('ModelConversionFcn',   defaultParams.ModelConversionFcn);
parser.addParameter('TransitionProbabilities', defaultParams.TransitionProbabilities);
parser.addParameter('ModelProbabilities',   defaultParams.ModelProbabilities);
parser.addParameter('MeasurementNoise',     defaultParams.MeasurementNoise);
parser.addParameter('ModelNames',           []);

% Accept heuristic NVs silently (already handled in constructor)
parser.addParameter('YawRateThresh',      []);
parser.addParameter('AccelThresh',        []);
parser.addParameter('Bias',               []);
parser.addParameter('Hysteresis',         []);
parser.addParameter('StickyDecay',        []);
parser.addParameter('UseWinnerTakesAll',  []);

% NEW optional NVs (if you want to expose them at construction):
parser.addParameter('RobustKind',         []); % 0/1/2
parser.addParameter('Nu',                 []);
parser.addParameter('NuJam',              []);
parser.addParameter('HuberK',             []);
parser.addParameter('HuberKJam',          []);
parser.addParameter('EpsSafeWTA',         []);
parser.addParameter('EpsInactiveMass',    []);
parser.addParameter('TopKMax',            []);

parser.parse(varargin{:});
args = parser.Results;

trackingFilters         = args.TrackingFilters;
transitionProbabilities = args.TransitionProbabilities;
modelConversion         = args.ModelConversionFcn;
modelProb               = args.ModelProbabilities;
measNoise               = args.MeasurementNoise;
state                   = args.State;
stateCov                = args.StateCovariance;
modelNames              = args.ModelNames;

% If user passed new NVs, pick them up (safe defaults already in object)
% You can set them post-construction as well.
end

% -------------------------------------------------------------------------
% Parse inputs for code generation
% -------------------------------------------------------------------------
function [state, stateCov, trackingFilters, transitionProbabilities, ...
          modelConversion, modelProb, measNoise, modelNames] = parseInputsCodegen(defaultParams, varargin)

parms = struct( ...
    'State',               uint32(0), ...
    'StateCovariance',     uint32(0), ...
    'TrackingFilters',     uint32(0), ...
    'ModelConversionFcn',  uint32(0), ...
    'TransitionProbabilities', uint32(0), ...
    'ModelProbabilities',  uint32(0), ...
    'MeasurementNoise',    uint32(0), ...
    'ModelNames',          uint32(0), ...
    'YawRateThresh',       uint32(0), ...
    'AccelThresh',         uint32(0), ...
    'Bias',                uint32(0), ...
    'Hysteresis',          uint32(0), ...
    'StickyDecay',         uint32(0), ...
    'UseWinnerTakesAll',   uint32(0), ...
    'RobustKind',          uint32(0), ...
    'Nu',                  uint32(0), ...
    'NuJam',               uint32(0), ...
    'HuberK',              uint32(0), ...
    'HuberKJam',           uint32(0), ...
    'EpsSafeWTA',          uint32(0), ...
    'EpsInactiveMass',     uint32(0), ...
    'TopKMax',             uint32(0));

popt = struct('CaseSensitivity', false, 'StructExpand', true, 'PartialMatching', false);

optarg = eml_parse_parameter_inputs(parms, popt, varargin{:});

trackingFilters         = eml_get_parameter_value(optarg.TrackingFilters,        defaultParams.TrackingFilters,        varargin{:});
transitionProbabilities = eml_get_parameter_value(optarg.TransitionProbabilities, defaultParams.TransitionProbabilities, varargin{:});
modelConversion         = eml_get_parameter_value(optarg.ModelConversionFcn,     defaultParams.ModelConversionFcn,     varargin{:});
modelProb               = eml_get_parameter_value(optarg.ModelProbabilities,     defaultParams.ModelProbabilities,     varargin{:});
measNoise               = eml_get_parameter_value(optarg.MeasurementNoise,       defaultParams.MeasurementNoise,       varargin{:});
state                   = eml_get_parameter_value(optarg.State,                  defaultParams.State,                  varargin{:});
stateCov                = eml_get_parameter_value(optarg.StateCovariance,        defaultParams.StateCovariance,        varargin{:});
modelNames              = coder.const(eml_get_parameter_value(optarg.ModelNames, [], varargin{:}));
end

