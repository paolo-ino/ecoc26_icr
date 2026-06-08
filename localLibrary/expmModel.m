classdef expmModel < MultiSectionalModel
    % Class for a multisectional model with sections having a transfer
    % matrix of the type expm(jB + C)

    properties
        beta0s
        dBeta0
        targetXtPerSection_dB
        variance_style % "uniform" or "exp_perMode": in the first case, the variance of the awgn components is the same across the whole matrix; in the second case it changes across groups. See details below.
        model_parameters % additional model parameters, like the coeff a used in exp(a*...) in the "exp_perMode" model
    end

    methods
        function obj = expmModel(beta0s, dBeta0, igsc, polFlag, groupSizes, modeIndices, nSec, varianceStyle, rndStream, varargin)
            obj@MultiSectionalModel(igsc, polFlag, groupSizes, modeIndices, nSec, rndStream);
            obj.beta0s = beta0s;
            obj.dBeta0 = dBeta0;
            obj.variance_style = varianceStyle;

            p = inputParser();
            p.addParameter("targetXtPerSection_dB", []) % it's good to have it empty since for functions like getTuningParameterFromXt. For the others, it being empty will result in an error, which is also the expected behavior
            p.addParameter("model_parameters", [])
            p.parse(varargin{:})

            obj.targetXtPerSection_dB = p.Results.targetXtPerSection_dB;

            obj.model_parameters = p.Results.model_parameters; % Empty for now. It is filled directly by generate_skew_hermitian_matrix_perMode if needed.
        end

        function H = singleSectionTranferMatrix(obj, modelParameters)
            targetXt_dB = obj.targetXtPerSection_dB;

            xt_dB_corrected = targetXt_dB+1-1.91+0.31*obj.nGroups; % correction to match the average XT (in the SingleValue metric sense)

            if obj.variance_style == "uniform"
                M = obj.generate_skew_hermitian_matrix();
                arg = 1j*diag(obj.beta0s) + sqrt( 1/(2*obj.nModes))*obj.dBeta0*sqrt(10.^(xt_dB_corrected/10))*M;
            elseif obj.variance_style == "exp_perMode"
                xt_dB_corrected = targetXt_dB+1-1.91+21+2.8; % correction to match the average XT (in the SingleValue metric sense)
                M = obj.generate_skew_hermitian_matrix_perMode(modelParameters{:});
                arg = 1j*diag(obj.beta0s) + sqrt( 1/(2*(obj.nModes)))*obj.dBeta0*sqrt(10.^(xt_dB_corrected/10))*M;            
            end

            H = expm(arg);
        end

        function H = generate_skew_hermitian_matrix_perMode(obj, varargin)
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
            
            B = 1/sqrt(2) * (randn(obj.rndStream, obj.nModes) + 1i*randn(obj.rndStream, obj.nModes));
        
            H = (B - B');
        end
    end

    methods (Static)
        function tuningParameter = getTuningParameterFromXt(obj, xt_dB)
            if obj.igsc == 1
                str = "_with_igsc";
            else
                str = "_no_igsc";
            end

            if obj.dBeta0 == 1
                str2 = "_dBeta0=1";
            elseif obj.dBeta0 == 1e4
                str2 = "_dBeta0=1e4";
            else
                str2 = "_dBeta0="+round(obj.dBeta0/1e4)+"e4";
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