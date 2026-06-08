classdef expmModel < MultiSectionalModel
%EXPMMODEL  Multi-sectional fibre model with matrix-exponential section transfer matrices.
%
%   expmModel implements a MULTISECTIONALMODEL where each section is described
%   by a transfer matrix of the form
%
%       H = expm(1j * diag(beta0sL) + C)
%
%   where diag(beta0sL) is the diagonal propagation-constant matrix multiplied 
%   by the section length (beta0sL can be regarded as an optimization parameter
%   or kept fixed at a certain reasonable value) and C is a
%   random skew-Hermitian coupling matrix. The variance of the elements of C
%   is scaled so that the average single-section crosstalk matches the target
%   TARGETXTPERSECTION_DB. Two variance styles are supported:
%
%     "uniform"     - the variance of every element of C is identical.
%     "exp_perMode" - the variance decays exponentially with the inter-group
%                     distance, producing stronger coupling between
%                     neighbouring mode groups.
%
%   Properties (own):
%     beta0sL              - [1 x nModes double] propagation constants of the
%                            fibre modes multiplied by the section length
%     dBeta0L              - [1 x 1 double] product between the difference in
%                            propagation constant among consecutive groups and
%                            the section length. It can be regarded as an
%                            optimisation parameter.
%     targetXtPerSection_dB- [1 x 1 double] target per-section crosstalk [dB];
%                            may be left empty for methods that do not require
%                            it (e.g. GETTUNINGPARAMETERFROMXT)
%     variance_style       - [string] "uniform" or "exp_perMode"; selects how
%                            the off-diagonal variance of C scales with mode
%                            group distance
%     model_parameters     - [struct] additional model parameters, e.g. the
%                            exponent a used in the "exp_perMode" scaling
%
%   Inherited properties: see MULTISECTIONALMODEL.
%
%   Example:
%     See example2_tuningParameter.m.
%
%   Author:      Paolo Carniello
%   Affiliation: Technical University of Munich
%   Date:        2026-06-08
%
%   See also MULTISECTIONALMODEL, XTHELPERS.

    properties
        beta0sL                % [1 x nModes double] propagation constants of the fibre modes (diagonal of B)  multiplied by the section length
        dBeta0L               % [1 x 1 double] product between the propagation-constant difference among consecutive groups and the section length
        targetXtPerSection_dB % [1 x 1 double] target per-section XT [dB]; can be empty when not needed
        variance_style        % [string] "uniform" or "exp_perMode": variance scaling of the coupling matrix elements
        model_parameters      % [struct] additional model parameters, e.g. the exponent a for the "exp_perMode" model
    end

    methods
        function obj = expmModel(beta0sL, dBeta0L, igsc, polFlag, groupSizes, modeIndices, nSec, varianceStyle, rndStream, varargin)
            %EXPMMODEL  Construct an expmModel.
            %
            %   OBJ = EXPMMODEL(BETA0SL, DBETA0L, IGSC, POLFLAG, GROUPSIZES,
            %       MODEINDICES, NSEC, VARIANCESTYLE, RNDSTREAM) creates an expmModel
            %       with the given parameters. TARGETXTPERSECTION_DB is left empty by
            %       default; pass it as a name-value argument when needed.
            %
            %   OBJ = EXPMMODEL(..., 'targetXtPerSection_dB', XT_DB) sets the target
            %       per-section crosstalk used in SINGLESECTIONTRANFERMATRIX.
            %
            %   OBJ = EXPMMODEL(..., 'model_parameters', PARAMS) sets extra model
            %       parameters (e.g. the exponent a for "exp_perMode").
            %
            %   Inputs:
            %     BETA0SL       - [1 x nModes double] propagation constants of the
            %                     fibre modes
            %     DBETA0L       - [1 x 1 double] product between the propagation-
            %                     constant difference among consecutive groups and the
            %                     section length
            %     IGSC          - [1 x 1 logical] enable intra-group strong coupling
            %     POLFLAG       - [1 x 1 logical] polarisation included in mode count
            %     GROUPSIZES    - [1 x nGroups double] number of modes per group
            %     MODEINDICES   - [1 x nGroups cell] mode index vectors per group
            %     NSEC          - [1 x 1 double] number of fibre sections
            %     VARIANCESTYLE - [string] "uniform" or "exp_perMode"
            %     RNDSTREAM     - [RandStream] random stream for all stochastic draws
            %
            %   Name-Value Arguments:
            %     'targetXtPerSection_dB' - [1 x 1 double] target per-section XT [dB];
            %                               default [] (empty)
            %     'model_parameters'      - [struct or double] extra model parameters;
            %                               default []
            %
            %   Output:
            %     OBJ - [expmModel] constructed object
            %
            %   Example:
            %     See example2_tuningParameter.m
            %
            %   See also MULTISECTIONALMODEL, SINGLESECTIONTRANFERMATRIX.
            obj@MultiSectionalModel(igsc, polFlag, groupSizes, modeIndices, nSec, rndStream);
            obj.beta0sL = beta0sL;
            obj.dBeta0L = dBeta0L;
            obj.variance_style = varianceStyle;

            p = inputParser();
            p.addParameter("targetXtPerSection_dB", []) % it's good to have it empty since for functions like getTuningParameterFromXt. For the others, it being empty will result in an error, which is also the expected behavior
            p.addParameter("model_parameters", [])
            p.parse(varargin{:})

            obj.targetXtPerSection_dB = p.Results.targetXtPerSection_dB;

            obj.model_parameters = p.Results.model_parameters; % Empty for now. It is filled directly by generate_skew_hermitian_matrix_perMode if needed.
        end

        function H = singleSectionTranferMatrix(obj, modelParameters)
            %SINGLESECTIONTRANFERMATRIX  Generate one section transfer matrix H = expm(1j*B + C).
            %
            %   H = SINGLESECTIONTRANFERMATRIX(OBJ, MODELPARAMETERS) returns a single
            %   realisation of the [NMODES x NMODES] transfer matrix
            %
            %       H = expm(1j * diag(beta0sL) + C)
            %
            %   where C is a random skew-Hermitian matrix scaled so that the average
            %   XT matches OBJ.TARGETXTPERSECTION_DB. The variance style
            %   OBJ.VARIANCE_STYLE selects whether C has uniform or inter-group-
            %   distance-dependent element variance.
            %
            %   Input:
            %     MODELPARAMETERS - [cell] extra arguments forwarded to
            %                       GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE;
            %                       used only for variance_style "exp_perMode"
            %
            %   Output:
            %     H - [nModes x nModes double] single-section transfer matrix
            %
            %   See also GENERATE_SKEW_HERMITIAN_MATRIX, GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE.
            targetXt_dB = obj.targetXtPerSection_dB;

            xt_dB_corrected = targetXt_dB-0.91+0.31*obj.nGroups; % correction to match the average XT (in the SingleValue metric sense)

            if obj.variance_style == "uniform"
                M = obj.generate_skew_hermitian_matrix();
                arg = 1j*diag(obj.beta0sL) + sqrt( 1/(2*obj.nModes))*obj.dBeta0L*sqrt(10.^(xt_dB_corrected/10))*M;
            elseif obj.variance_style == "exp_perMode"
                xt_dB_corrected = targetXt_dB+22.89; % correction to match the average XT (in the SingleValue metric sense)
                M = obj.generate_skew_hermitian_matrix_perMode(modelParameters{:});
                arg = 1j*diag(obj.beta0sL) + sqrt( 1/(2*(obj.nModes)))*obj.dBeta0L*sqrt(10.^(xt_dB_corrected/10))*M;
            end

            H = expm(arg);
        end

        function H = generate_skew_hermitian_matrix_perMode(obj, varargin)
            %GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE  Random skew-Hermitian matrix with group-distance-dependent variance.
            %
            %   H = GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE(OBJ) generates an
            %   [NMODES x NMODES] skew-Hermitian matrix whose element variance decays
            %   exponentially with the inter-group distance:
            %
            %       stddev(H(i,j)) ~ exp(|groupIndex(i) - groupIndex(j)|)^a
            %
            %   This models stronger coupling between neighbouring mode groups.
            %
            %   H = GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE(OBJ, 'a', A) uses the
            %   exponent A (default -2.25).
            %
            %   Name-Value Arguments:
            %     'a' - [1 x 1 double] exponent for the inter-group variance decay;
            %           default -2.25
            %
            %   Output:
            %     H - [nModes x nModes double] skew-Hermitian coupling matrix
            %
            %   See also GENERATE_SKEW_HERMITIAN_MATRIX, SINGLESECTIONTRANFERMATRIX.
            nGroups = length(obj.modeIndices);
            nModes = 0;
            for gg = 1:nGroups
                nModes = nModes + length(obj.modeIndices{gg});
            end

            varianceMat = zeros(nModes, nModes);
            for rowGroupIdx = 1:nGroups
                for colGroupIdx = 1:nGroups
                    switch obj.variance_style
                        case "exp_perMode"
                            p = inputParser();
                            p.addParameter("a", -2.25)
                            p.parse(varargin{:})

                            a = p.Results.a;
                            obj.model_parameters.a = a;
                            scaling = exp(abs(colGroupIdx - rowGroupIdx))^a;
                            matBlock = scaling.^2 * ones(length(obj.modeIndices{rowGroupIdx}), length(obj.modeIndices{colGroupIdx}));
                        otherwise
                            error("Unsupported variance scaling")
                    end
                    varianceMat(obj.modeIndices{rowGroupIdx}, obj.modeIndices{colGroupIdx}) = matBlock;
                end
            end

            B = 1/sqrt(2) * (randn(obj.rndStream, nModes) + 1i*randn(obj.rndStream, nModes));
            B = B .* sqrt(varianceMat); % scale the variance of the various elements according to the variance matrix
            H = (B - B');
        end

        function H = generate_skew_hermitian_matrix(obj)
            %GENERATE_SKEW_HERMITIAN_MATRIX  Random skew-Hermitian matrix with uniform variance.
            %
            %   H = GENERATE_SKEW_HERMITIAN_MATRIX(OBJ) generates an
            %   [NMODES x NMODES] skew-Hermitian matrix with identical variance for
            %   every element (the "uniform" variance style).
            %
            %   Output:
            %     H - [nModes x nModes double] skew-Hermitian coupling matrix
            %
            %   See also GENERATE_SKEW_HERMITIAN_MATRIX_PERMODE, SINGLESECTIONTRANFERMATRIX.

            B = 1/sqrt(2) * (randn(obj.rndStream, obj.nModes) + 1i*randn(obj.rndStream, obj.nModes));

            H = (B - B');
        end
    end

    methods (Static)
        function tuningParameter = getTuningParameterFromXt(obj, xt_dB)
            %GETTUNINGPARAMETERFROMXT  Map a target XT [dB] to the corresponding tuning parameter.
            %
            %   TUNINGPARAMETER = GETTUNINGPARAMETERFROMXT(OBJ, XT_DB) loads the
            %   precomputed interpolant matching the model configuration (igsc,
            %   dBeta0L, polFlag, variance_style) from the databases folder and
            %   evaluates it at XT_DB.
            %
            %   The database file is selected automatically from the properties of OBJ.
            %   Use MULTISECTIONALMODEL.COMPUTEXTVSSTUNINGPARAMETER and
            %   MULTISECTIONALMODEL.BUILDINTERPOLANT to generate the database files.
            %
            %   Inputs:
            %     OBJ    - [expmModel] model object whose properties select the
            %              database file (igsc, dBeta0L, polFlag, variance_style,
            %              model_parameters, nGroups)
            %     XT_DB  - [1 x 1 double] target crosstalk [dB]
            %
            %   Output:
            %     TUNINGPARAMETER - [1 x 1 double] tuning-parameter value
            %                       (i.e. targetXtPerSection_dB) that achieves XT_DB
            %
            %   See also MULTISECTIONALMODEL.BUILDINTERPOLANT,
            %            MULTISECTIONALMODEL.COMPUTEXTVSSTUNINGPARAMETER.
            if obj.igsc == 1
                str = "_with_igsc";
            else
                str = "_no_igsc";
            end

            if obj.dBeta0L == 1
                str2 = "_dBeta0=1";
            elseif obj.dBeta0L == 1e4
                str2 = "_dBeta0=1e4";
            else
                str2 = "_dBeta0="+round(obj.dBeta0L/1e4)+"e4";
            end

            if obj.polFlag == 1
                str3 = "_with_pol";
            else
                str3 = "_no_pol";
            end

            if obj.variance_style == "uniform"
                str4 = "_unif";
            elseif obj.variance_style == "exp_perMode"
                str4 = "_exp_per_mode_a="+obj.model_parameters.a;
            end

            load("databases/single_segment_dB_vs_tuning_parameter_expm"+str4+str2+str3+str, "interpolants") % all parameters: "nGroupsVec", "tuningParamName", "tuningParamVec", "xtSingleValueVec", "sectionObj", "xtMetricsArray", "interpolants"

            tuningParameter = interpolants{obj.nGroups}(xt_dB);
        end
    end
end
