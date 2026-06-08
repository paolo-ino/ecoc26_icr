classdef (Abstract) MultiSectionalModel < handle
%MULTISECTIONALMODEL  Abstract base class for multi-sectional optical fibre models.
%
%   MultiSectionalModel provides the common framework for simulating mode
%   coupling in few-mode or multimode optical fibres using a transfer-matrix
%   approach. A fibre is modelled as the concatenation of NSEC independent
%   sections, each described by an [NMODES x NMODES] transfer matrix.
%   Concrete subclasses (e.g. EXPMMODEL) must implement
%   SINGLESECTIONTRANFERMATRIX, which generates one realisation of a
%   single-section transfer matrix, and GETTUNINGPARAMETERFROMXT, which
%   inverts the XT-vs-parameter relation.
%
%   The class provides:
%     - Monte-Carlo section concatenation    (CONCATENATIONTRANSFERMATRIX)
%     - Single-section statistics            (SINGLESECTIONSTATISTICS)
%     - Intra-group strong-coupling mixing   (INTRAGROUPSTRONGCOUPLING)
%     - Crosstalk computation                (COMPUTEXT)
%     - XT-vs-tuning-parameter sweep (static)(COMPUTEXTVSSTUNINGPARAMETER)
%     - Interpolant construction     (static)(BUILDINTERPOLANT)
%
%   MultiSectionalModel inherits from handle so that property values persist
%   across methods on the same instance.
%
%   Note: in subclass constructors, pass the tuning parameter as an optional
%   name-value argument so it can be swept via COMPUTEXTVSSTUNINGPARAMETER.
%
%   Properties:
%     igsc        - [1 x 1 logical] flag enabling intra-group strong coupling
%                   through Haar-random unitary matrices applied after each
%                   section
%     groupSizes  - [1 x nGroups double] number of modes in each mode group
%     modeIndices - [1 x nGroups cell] each cell element is a
%                   [1 x nModesPerGroup double] vector of global mode indices
%                   belonging to that group
%     nGroups     - [1 x 1 double] number of mode groups; derived
%                   automatically from groupSizes
%     nModes      - [1 x 1 double] total number of modes; derived
%                   automatically as sum(groupSizes)
%     polFlag     - [1 x 1 logical] flag indicating that polarisation degrees
%                   of freedom are included in the mode count
%     nSec        - [1 x 1 double] number of fibre sections in the
%                   concatenation
%     rndStream   - [RandStream] random stream controlling all stochastic
%                   draws; create as, e.g., RandStream('mrg32k3a', 'Seed', yourSeed)
%
%   Example:
%     % See example2_concatenation.m for a concrete instantiation example.
%
%   Author:      Paolo Carniello
%   Affiliation: Technical University of Munich
%   Date:        2026-06-08
%
%   See also EXPMMODEL, XTHELPERS.

    properties
        igsc        % [1 x 1 logical] flag for intra-group strong coupling through random unitary matrices
        groupSizes  % [1 x nGroups double] number of modes per group
        modeIndices % [1 x nGroups cell] each element is a [1 x nModesPerGroup double] array of mode indices per group
        nGroups     % [1 x 1 double] number of mode groups (derived from groupSizes)
        nModes      % [1 x 1 double] total number of modes (derived as sum(groupSizes))
        polFlag     % [1 x 1 logical] flag indicating polarisation is included
        nSec        % [1 x 1 double] number of fibre sections in the concatenation
        rndStream   % [RandStream] random stream; create as RandStream('mrg32k3a', 'Seed', seed)
    end

    methods
        function obj = MultiSectionalModel(igsc, polFlag, groupSizes, modeIndices, nSec, rndStream)
            %MULTISECTIONALMODEL  Construct a MultiSectionalModel and initialise shared parameters.
            %
            %   OBJ = MULTISECTIONALMODEL(IGSC, POLFLAG, GROUPSIZES, MODEINDICES,
            %       NSEC, RNDSTREAM) stores the model parameters common to all
            %       concrete subclasses. NGROUPS and NMODES are derived automatically
            %       from GROUPSIZES.
            %
            %   Inputs:
            %     IGSC        - [1 x 1 logical] enable intra-group strong coupling
            %     POLFLAG     - [1 x 1 logical] polarisation included in mode count
            %     GROUPSIZES  - [1 x nGroups double] number of modes per group
            %     MODEINDICES - [1 x nGroups cell] mode index vectors per group
            %     NSEC        - [1 x 1 double] number of fibre sections
            %     RNDSTREAM   - [RandStream] random stream for all stochastic draws
            %
            %   Output:
            %     OBJ - [MultiSectionalModel] constructed object handle
            %
            %   See also EXPMMODEL.
            obj.igsc = igsc;
            obj.groupSizes = groupSizes;
            obj.modeIndices = modeIndices;
            obj.nSec = nSec;
            obj.rndStream = rndStream;
            obj.nGroups = length(obj.groupSizes);
            obj.nModes = sum(obj.groupSizes);
            obj.polFlag = polFlag;
        end

        function [T_tot, xtMetrics] = concatenationTransferMatrix(obj, varargin)
            %CONCATENATIONTRANSFERMATRIX  Monte-Carlo concatenation of fibre sections.
            %
            %   [T_TOT, XTMETRICS] = CONCATENATIONTRANSFERMATRIX(OBJ, 'nRep', NREP)
            %   runs NREP Monte-Carlo repetitions of a multi-sectional fibre. In each
            %   repetition, NSEC section transfer matrices are generated via
            %   SINGLESECTIONTRANFERMATRIX, optionally post-multiplied by an
            %   intra-group strong-coupling matrix, and cascaded sequentially.
            %   Crosstalk metrics are evaluated on the cumulative transfer matrix
            %   after each section.
            %
            %   [T_TOT, XTMETRICS] = CONCATENATIONTRANSFERMATRIX(OBJ, ..., Name, Value)
            %   accepts the following optional name-value arguments:
            %
            %   Name-Value Arguments:
            %     'xtMetrics'        - [1 x nXtDefs cell of strings] XT definitions
            %                          to compute; default {"Ferreira"}
            %     'nRep'             - [1 x 1 double] number of Monte-Carlo repetitions
            %     'singleSectionArgs'- [cell] extra arguments forwarded to
            %                          SINGLESECTIONTRANFERMATRIX; default {}
            %
            %   Outputs:
            %     T_TOT     - [nRep x nSec x nModes x nModes double] cumulative
            %                 transfer matrices; T_TOT(r,s,:,:) is the product of
            %                 sections 1..s in repetition r
            %     XTMETRICS - [1 x nXtDefs struct] array with fields:
            %                   perGroup - [nRep x nSec x nGroups] per-group XT
            %                   avgDef1  - [nRep x nSec] primary averaged XT
            %                   avgDef2  - [nRep x nSec] secondary averaged XT
            %                             (NaN if unused by the definition)
            %
            %   Example:
            %     [T, xt] = obj.concatenationTransferMatrix('nRep', 100, ...
            %                   'xtMetrics', {"Ferreira", "SingleValue"});
            %
            %   See also SINGLESECTIONSTATISTICS, COMPUTEXT.
            p = inputParser();
            p.addParameter("xtMetrics", {"Ferreira"}) % cell array of strings with the names of the xt definitions to use for the computation of xt.
            p.addParameter("nRep", []) % nr of monte-carlo repetitions of the concatenation
            p.addParameter("singleSectionArgs", {}) % the function can be called without passing extra args for the single section function
            p.parse(varargin{:})

            nXtDefs = length(p.Results.xtMetrics);
            nRep = p.Results.nRep;

            T_tot = zeros(nRep, obj.nSec, obj.nModes, obj.nModes);
            xtMetrics = struct();

            for repIdx = 1:nRep
                disp("Rep: "+repIdx + "/" + nRep)
                H_overall = eye(obj.nModes);
                for secIdx = 1:obj.nSec
                    %% Generate a single segment
                    H_segment = obj.singleSectionTranferMatrix(p.Results.singleSectionArgs);

                    if obj.igsc
                        % Force strong intra-group coupling
                        R = obj.intraGroupStrongCoupling();
                        H_segment = H_segment*R;
                    end

                    % Compute the overall transfer matrix
                    H_overall = H_segment * H_overall;
                    T_tot(repIdx, secIdx, :, :) = H_overall;

                    P_mat_dBm = 10*log10(abs(H_overall).^2);

                    %% Compute XT
                    for xtDefIdx = 1:nXtDefs
                        [outPerGroup, avg1, avg2] =...
                            obj.computeXT(P_mat_dBm, p.Results.xtMetrics{xtDefIdx});

                        % Depending on the xt metric, some of the outputs are
                        % empty. To be able to store them in an array, you need
                        % to do something about it, e.g., save them as NaNs.
                        if isempty(outPerGroup); outPerGroup = NaN; end
                        if isempty(avg1); avg1 = NaN; end
                        if isempty(avg2); avg2 = NaN; end

                        xtMetrics(xtDefIdx).perGroup(repIdx, secIdx, :) = outPerGroup;
                        xtMetrics(xtDefIdx).avgDef1(repIdx, secIdx) = avg1;
                        xtMetrics(xtDefIdx).avgDef2(repIdx, secIdx) = avg2;
                    end
                end
            end
        end

        function R = intraGroupStrongCoupling(obj)
            %INTRAGROUPSTRONGCOUPLING  Block-diagonal Haar-random unitary mixing matrix.
            %
            %   R = INTRAGROUPSTRONGCOUPLING(OBJ) builds an [NMODES x NMODES]
            %   block-diagonal unitary matrix. Each diagonal block is an independent
            %   Haar-distributed random unitary matrix of size equal to the
            %   corresponding group, enforcing intra-group mode mixing while
            %   leaving inter-group coupling unaffected.
            %
            %   Output:
            %     R - [nModes x nModes double] block-diagonal unitary matrix
            %         (one Haar-random unitary block per mode group)
            %
            %   See also XTHELPERS.RANDOMUNIATRYMATRIX, SINGLESECTIONSTATISTICS.
            R1 = cell(1, length(obj.groupSizes));
            for groupIdx = 1:length(obj.groupSizes)
                R1{groupIdx} = xtHelpers.randomUnitaryMatrix(obj.groupSizes(groupIdx), obj.rndStream); % polarization coupling
            end
            R = blkdiag(R1{:});
        end

        function [xtInterPerGroup, xtInterAvgDef1, xtInterAvgDef2] = computeXT(obj, P_mat_dB, def)
            %COMPUTEXT  Compute inter-group crosstalk from a power-coupling matrix.
            %
            %   [XTINTERPERGROUP, XTINTERAVGDEF1, XTINTERAVGDEF2] = ...
            %       COMPUTEXT(OBJ, P_MAT_DB, DEF) computes inter-group XT according
            %       to the definition DEF.
            %
            %   Supported definitions:
            %     "Ferreira"    - per-group XT from [Ferreira et al., 2017]
            %     "SingleValue" - single scalar XT; reciprocal of the metric in
            %                     [Mazur et al., 2019]
            %
            %   Inputs:
            %     P_MAT_DB - [nModes x nModes double] power-coupling matrix in dB
            %     DEF      - [string] XT definition: "Ferreira" or "SingleValue"
            %
            %   Outputs:
            %     XTINTERPERGROUP - [1 x nGroups double] per-group XT in linear scale;
            %                       [] for "SingleValue"
            %     XTINTERAVGDEF1  - [1 x 1 double] primary aggregated XT (linear)
            %     XTINTERAVGDEF2  - [1 x 1 double] secondary aggregated XT (linear);
            %                       [] for the currently supported definitions
            %
            %   References:
            %     [1] Ferreira et al., "Semi-Analytical Modelling of Linear Mode
            %         Coupling in Few-Mode Fibers", 2017.
            %     [2] Mazur et al., "Characterization of Long Multi-Mode Fiber Links
            %         using Digital Holography", 2019.
            %
            %   See also SINGLESECTIONSTATISTICS, CONCATENATIONTRANSFERMATRIX.
            P_mat = 10.^(P_mat_dB/10);

            nModesPerGroup = zeros(1, obj.nGroups);
            for g = 1:obj.nGroups
                nModesPerGroup(g) = length(obj.modeIndices{g});
            end

            xtInterPerGroup = zeros(1, obj.nGroups);
            switch def
                case "Ferreira"
                    % From [Ferreira et al., "Semi-Analytical Modelling of Linear
                    % Mode Coupling in Few-Mode Fibers", 2017]
                    for gIdx = 1:obj.nGroups
                        Psub = P_mat(:, obj.modeIndices{gIdx}); % vertical stripe
                        Psub_sum = sum(sum(Psub));
                        Pblock = sum(sum(P_mat(obj.modeIndices{gIdx}, obj.modeIndices{gIdx})));
                        xtInterPerGroup(gIdx) = (Psub_sum-Pblock)./Pblock;
                    end
                    xtInterAvgDef1 = sum(xtInterPerGroup.*nModesPerGroup)/obj.nModes;
                    xtInterAvgDef2 = [];
                case "SingleValue"
                    % Reciprocal of the quantity defined in [Mazur et al.,
                    % "Characterization of Long Multi-Mode Fiber Links using Digital Holography", 2019]
                    Pblock = zeros(1, obj.nGroups);
                    for gIdx = 1:obj.nGroups
                        Pblock(gIdx) = sum(sum(P_mat(obj.modeIndices{gIdx}, obj.modeIndices{gIdx})));
                    end
                    xtInterSingleValue = sum(Pblock)./(sum(sum(P_mat))-sum(Pblock)); % []
                    xtInterAvgDef1 = 1/xtInterSingleValue;
                    xtInterPerGroup = [];
                    xtInterAvgDef2 = [];
            end
        end

        function [H_vec, xtMetrics] = singleSectionStatistics(obj, nRep, varargin)
            %SINGLESECTIONSTATISTICS  Monte-Carlo statistics of a single fibre section.
            %
            %   [H_VEC, XTMETRICS] = SINGLESECTIONSTATISTICS(OBJ, NREP) generates NREP
            %   independent realisations of a single-section transfer matrix via
            %   SINGLESECTIONTRANFERMATRIX, optionally post-multiplied by an
            %   intra-group strong-coupling matrix, and returns the realisations
            %   together with their crosstalk metrics.
            %
            %   [H_VEC, XTMETRICS] = SINGLESECTIONSTATISTICS(OBJ, NREP, Name, Value)
            %   accepts the following optional name-value arguments:
            %
            %   Name-Value Arguments:
            %     'xtMetrics'        - [1 x nXtDefs cell of strings] XT definitions
            %                          to compute; default {"Ferreira"}
            %     'singleSectionArgs'- [cell] extra arguments forwarded to
            %                          SINGLESECTIONTRANFERMATRIX; default {}
            %
            %   Inputs:
            %     NREP - [1 x 1 double] number of Monte-Carlo realisations
            %
            %   Outputs:
            %     H_VEC     - [nRep x nModes x nModes double] single-section
            %                 transfer matrices
            %     XTMETRICS - [1 x nXtDefs struct] array with fields:
            %                   perGroup - [nRep x nGroups] per-group XT
            %                   avgDef1  - [nRep x 1] primary averaged XT
            %                   avgDef2  - [nRep x 1] secondary averaged XT
            %                             (NaN if unused by the definition)
            %
            %   Example:
            %     [H, xt] = obj.singleSectionStatistics(100, ...
            %                   'xtMetrics', {"Ferreira", "SingleValue"});
            %
            %   See also CONCATENATIONTRANSFERMATRIX, COMPUTEXT.
            p = inputParser();
            p.addParameter("xtMetrics", {"Ferreira"}) % cell array of strings with the names of the xt definitions to use for the computation of xt.
            p.addParameter("singleSectionArgs", {}) % parameters to pass to singleSectionTranferMatrix
            p.parse(varargin{:})

            nXtDefs = length(p.Results.xtMetrics);

            H_vec = zeros(nRep, obj.nModes, obj.nModes);
            % xtMetrics = struct('perGroup', zeros(nRep, nGroups), 'avgDef1', zeros(1, nRep), 'avgDef2', zeros(1, nRep));
            xtMetrics = struct();

            for repIdx = 1:nRep
                %% Generate a single segment
                H = obj.singleSectionTranferMatrix(p.Results.singleSectionArgs);

                if obj.igsc
                    % Force strong intra-group coupling
                    R = obj.intraGroupStrongCoupling();
                    H = H*R;
                end
                P_mat_dBm = 10*log10(abs(H).^2);

                %% Compute XT

                for xtDefIdx = 1:nXtDefs
                    [outPerGroup, avg1, avg2] = obj.computeXT(P_mat_dBm, p.Results.xtMetrics{xtDefIdx});

                    % Depending on the xt metric, some of the outputs are
                    % empty. To be able to store them in an array, you need
                    % to do something about it, e.g., save them as NaNs.
                    if isempty(outPerGroup); outPerGroup = NaN; end
                    if isempty(avg1); avg1 = NaN; end
                    if isempty(avg2); avg2 = NaN; end

                    xtMetrics(xtDefIdx).perGroup(repIdx, :) = outPerGroup;
                    xtMetrics(xtDefIdx).avgDef1(repIdx) = avg1;
                    xtMetrics(xtDefIdx).avgDef2(repIdx) = avg2;
                end

                H_vec(repIdx, :, :) = H;
            end
        end
    end

    methods (Static)
        function [xtAvgSingleValueLoop, xtMetrics, lastObj] = computeXtVsTuningParameter(ModelConstructor, tuningParamName, tuningParamVec, ModelConstructorArgsWihoutTuningParam, xtMetricsArray, nRep, varargin)
            %COMPUTEXTVSSTUNINGPARAMETER  Sweep a model tuning parameter and record the resulting XT.
            %
            %   [XTAVGSINGLEVALUELOOP, XTMETRICS, LASTOBJ] = ...
            %       COMPUTEXTVSSTUNINGPARAMETER(MODELCONSTRUCTOR, TUNINGPARAMNAME, ...
            %       TUNINGPARAMVEC, MODELCONSTRUCTORARGSWITHOUTTUNINGPARAM, ...
            %       XTMETRICSARRAY, NREP) sweeps TUNINGPARAMVEC, builds a model at
            %       each value via MODELCONSTRUCTOR, and calls SINGLESECTIONSTATISTICS
            %       with NREP realisations to estimate the resulting XT.
            %
            %   [...] = COMPUTEXTVSSTUNINGPARAMETER(..., 'singleSectionArgs', ARGS)
            %   forwards ARGS to SINGLESECTIONSTATISTICS.
            %
            %   Inputs:
            %     MODELCONSTRUCTOR                       - [function handle] constructor
            %                                              of the concrete model class,
            %                                              e.g. @expmModel
            %     TUNINGPARAMNAME                        - [string] name-value key for
            %                                              the tuning parameter passed to
            %                                              the constructor, e.g.
            %                                              "targetXtPerSection_dB"
            %     TUNINGPARAMVEC                         - [1 x N double] values of the
            %                                              tuning parameter to sweep
            %     MODELCONSTRUCTORARGSWITHOUTTUNINGPARAM - [cell] all constructor
            %                                              arguments excluding the tuning
            %                                              parameter, e.g.
            %                                              {beta0s, dBeta0L, igsc, ...}
            %     XTMETRICSARRAY                         - [1 x nXtDefs cell of strings]
            %                                              XT definitions to compute
            %     NREP                                   - [1 x 1 double] number of
            %                                              single-section realisations
            %                                              per tuning value
            %
            %   Name-Value Arguments:
            %     'singleSectionArgs' - [cell] extra arguments forwarded to
            %                           SINGLESECTIONSTATISTICS; default {}
            %
            %   Outputs:
            %     XTAVGSINGLEVALUELOOP - [1 x N double] mean "SingleValue" XT for each
            %                            element of TUNINGPARAMVEC (linear scale)
            %     XTMETRICS            - [nXtDefs x N struct] raw per-realisation
            %                            metrics with fields perGroup, avgDef1, avgDef2
            %     LASTOBJ              - [MultiSectionalModel] model object built for
            %                            the last tuning-parameter value; returned for
            %                            convenient parameter inspection and saving
            %
            %   Example:
            %     See example3_tuningParameter.m for a concrete example.
            %
            %   See also BUILDINTERPOLANT, SINGLESECTIONSTATISTICS.
            p = inputParser();
            p.addParameter("singleSectionArgs", {}) % extra optional arguments to pass to singleSectionStatistics
            p.parse(varargin{:})

            singleSectionArgs = p.Results.singleSectionArgs;

            paramLoopLength = numel(tuningParamVec);
            nXtDef = numel(xtMetricsArray);

            for tuningParamIdx = 1:paramLoopLength
                args = [ModelConstructorArgsWihoutTuningParam, {tuningParamName, tuningParamVec(tuningParamIdx)}];
                lastObj = ModelConstructor(args{:});

                disp("Tuning param loop: "+tuningParamIdx +"/"+paramLoopLength)
                % You cannot preallocate them differently, as the fields have
                % different sizes as nGroups changes
                xtMetrics(nXtDef, tuningParamIdx) = struct('perGroup', zeros(nRep, lastObj.nGroups), 'avgDef1', zeros(1, nRep), 'avgDef2', zeros(1, nRep));
                xtMetricsAvg(nXtDef, tuningParamIdx) = struct('perGroup', zeros(1, lastObj.nGroups), 'avgDef1', zeros(1, 1), 'avgDef2', zeros(1, 1));
                [~, xtMetrics(:, tuningParamIdx)] = lastObj.singleSectionStatistics(nRep, "xtMetrics", xtMetricsArray, "singleSectionArgs", singleSectionArgs);

                for xtDefIdx = 1:nXtDef
                    xtMetricsAvg(xtDefIdx, tuningParamIdx).perGroup = mean(xtMetrics(xtDefIdx, tuningParamIdx).perGroup, 1);
                    xtMetricsAvg(xtDefIdx, tuningParamIdx).avgDef1 =  mean(xtMetrics(xtDefIdx, tuningParamIdx).avgDef1);
                    xtMetricsAvg(xtDefIdx, tuningParamIdx).avgDef2 =  mean(xtMetrics(xtDefIdx, tuningParamIdx).avgDef2);
                end
            end

            xtAvgSingleValueLoop = zeros(1, paramLoopLength);
            for tuningParamIdx = 1:paramLoopLength
                [~, loc] = ismember("SingleValue", string(xtMetricsArray));
                xtAvgSingleValueLoop(tuningParamIdx) = xtMetricsAvg(loc, tuningParamIdx).avgDef1;
            end
        end

        function interpolant = buildInterpolant(tuningParamVec, xtVec, windowSize)
            %BUILDINTERPOLANT  Build a monotonic interpolant mapping XT [dB] to a tuning parameter.
            %
            %   INTERPOLANT = BUILDINTERPOLANT(TUNINGPARAMVEC, XTVEC, WINDOWSIZE)
            %   constructs a function handle that maps a target XT value in dB to the
            %   corresponding tuning-parameter value.
            %
            %   The XT curve is first converted to dB. A naive moving average of window
            %   WINDOWSIZE is applied to the tail of the curve to enforce strict
            %   monotonicity (required for bijectivity). A linear griddedInterpolant
            %   is then fitted on the monotonic region. The returned function applies
            %   linear extrapolation for inputs below the training range and saturates
            %   at the maximum tuning-parameter value for inputs above it.
            %
            %   Inputs:
            %     TUNINGPARAMVEC - [1 x N double] tuning-parameter values
            %     XTVEC          - [1 x N double] corresponding XT values in linear
            %                      scale (not dB); must be positive
            %     WINDOWSIZE     - [1 x 1 double] moving-average window length used
            %                      to enforce monotonicity in the XT-dB curve
            %
            %   Output:
            %     INTERPOLANT - [function handle] maps xt_dB_input ([1 x M double])
            %                   to the corresponding tuning-parameter values
            %                   ([1 x M double])
            %
            %   Example:
            %     F = MultiSectionalModel.buildInterpolant(tuningVec, xtVec, 5);
            %     targetParam = F(-20);   % tuning parameter yielding XT = -20 dB
            %
            %   See also COMPUTEXTVSSTUNINGPARAMETER, GRIDDEDINTERPOLANT.
            xtVec_dB = 10*log10(xtVec);

            % {
            maxIdx = find(diff(xtVec_dB)<=0, 1, 'first'); % to ensure the mapping to be bijective

            if isempty(maxIdx)
                maxIdx = size(xtVec_dB, 2);
            else
                ff = movmean(xtVec_dB(maxIdx-windowSize+1:end), windowSize);

                xtVec_dB(maxIdx+1-windowSize:end) = ff;
                maxIdx = find(diff(xtVec_dB)<=0, 1, 'first'); % to ensure the mapping to be bijective
                if isempty(maxIdx)
                    maxIdx = size(xtVec_dB, 2);
                end
            end
            %}

            F = griddedInterpolant(xtVec_dB(1:maxIdx), tuningParamVec(1:maxIdx), 'linear', 'linear');  % linear inside + linear extrap

            xmin = min(xtVec_dB(1:maxIdx));
            xmax = max(xtVec_dB(1:maxIdx));
            ymax = tuningParamVec(maxIdx);   % assuming tuningParamVec is sorted

            % Ensure that for inputs which are below the training set, we follow a
            % linear behavior (cause that's what I expect/know to be correct), and
            % above apply a saturation.
            interpolant = @(xt_dB_input) arrayfun(@(xi) ...
                (xi < xmin) * F(xi) + ...
                (xi >= xmin && xi <= xmax) * F(xi) + ...
                (xi > xmax) * ymax, ...
                xt_dB_input);
        end
    end

    methods (Abstract)
        %SINGLESECTIONTRANFERMATRIX  Generate one [nModes x nModes] section transfer matrix.
        %   Must be implemented by each concrete subclass.
        H = singleSectionTranferMatrix(obj, varargin);

        %GETTUNINGPARAMETERFROMXT  Map a target XT [dB] to the tuning-parameter value.
        %   Must be implemented by each concrete subclass.
        tuningParameter = getTuningParameterFromXt(obj, xt_dB, varargin);
    end
end
