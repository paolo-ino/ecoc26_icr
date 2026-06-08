classdef (Abstract) MultiSectionalModel < handle
    % MultiSectionalModel is of type handle to be able to save some
    % parameters across methods.
    % 
    % Some hints:
    % - when creating a child class, make sure that the tuning parameter is
    % passed as optional argument to the constructor, so that you can loop
    % over that
    properties
        igsc % flag for intra-group strong coupling through random unitary matrices
        groupSizes
        modeIndices  % [1 x nGroups] cell array, where each element is an [1 x nModesPerGroup] array with the indices of the modes which belong to the group
        nGroups
        nModes
        polFlag
        nSec
        rndStream % Define rndStream similar to rndStream = RandStream('mrg32k3a', 'Seed', 2007509295)
    end

    methods
        function obj = MultiSectionalModel(igsc, polFlag, groupSizes, modeIndices, nSec, rndStream)
            % Not sure what to do here
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
            p = inputParser();
            p.addParameter("xtMetrics", {"Ferreira"}) % cell array of strings with the names of the xt definitions to use for the computation of xt.
            p.addParameter("nRep", []) % nr of monte-carlo repetitions of the concatenation
            p.addParameter("singleSectionArgs", {}) % the function can be called without passing extra args for the single section function
            p.parse(varargin{:})

            nXtDefs = length(p.Results.xtMetrics);
            nRep = p.Results.nRep;

            T_tot = zeros(nRep, obj.nSec, obj.nModes, obj.nModes);
            % xtMetrics = struct('perGroup', zeros(nRep, nGroups), 'avgDef1', zeros(1, nRep), 'avgDef2', zeros(1, nRep));
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
            % Create a block-diagonal matrix R whose blocks are random
            % unitary matrices uniformly distributed over U(N) according to
            % the Haar metric.
            R1 = cell(1, length(obj.groupSizes));
            for groupIdx = 1:length(obj.groupSizes)
                R1{groupIdx} = xtHelpers.randomUnitaryMatrix(obj.groupSizes(groupIdx), obj.rndStream); % polarization coupling
            end
            R = blkdiag(R1{:});
        end

        function [xtInterPerGroup, xtInterAvgDef1, xtInterAvgDef2] = computeXT(obj, P_mat_dB, def)
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

        % A method which produces nRep realizations of the TM for a
        % section. It accepts a callback function to compute a single
        % realization of the TM of a section. It allows to include
        % intra-group strong coupling through random unitary matrices. It
        % accepts a cell array of strings with the names of the xt
        % definitions to use for the computation of xt.
        function [H_vec, xtMetrics] = singleSectionStatistics(obj, nRep, varargin)

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
            % Compute the relation between XT (according to the
            % "SingleValue" metric) and the tuning parameter of the model.
            %
            % INPUTS:
            % ModelConstructor: eg: @expmModel
            % ModelConstructorArgs: eg: {nModes, ...
            %                           nGroups, ...
            %                           beta0s, ...
            %                           dBeta0, ...
            %                           targetXtPerSegment_dB_vec(xtIdx), ...
            %                           rndStream, ...
            %                           "varianceStyle", varianceStyle, ...
            %                           "modeIndices", modeIndices};
            %
            % OUTPUTS:
            % lastObj: returning the last object so one can save more
            % easily the object parameters.

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
            % Build the interpolation relation between XT and the tuning
            % parameter of the model.
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
            ymax = tuningParamVec(maxIdx);   % assuming x is sorted

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
        H = singleSectionTranferMatrix(obj, varargin);

        tuningParameter = getTuningParameterFromXt(obj, xt_dB, varargin);
    end
end