% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Example script showing how to compute the relation between the tuning
% parameter of the model and the target XT.
% 
% The script has already been run for the cases presented in the paper. The
% relations have been saved in the folder ./databases/. The script can be
% rerun for additional cases.
% To see how to use the numerical relations, check
% example1_single_section.m or example2_concatenation.m.
% 
% 
% References:
% 
% [Carniello26] Carniello et al., "A Simplified Model for Linear Mode Coupling in
% Multimode Fibers with Experimental Assessment", submitted to ECOC 2026
% 
%   Author:      Paolo Carniello
%   Affiliation: Technical University of Munich
%   Date:        2026-06-08
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clearvars

addpath(genpath("./localLibrary"))
addpath(genpath("./databases"))

%% Parameters you might want to change
% nGroupsVec = [3 9];
nGroupsVec = 2:12;
polFlag = 1; % flag to decide if polarizations are considered as well, or only spatial modes
nRep = 100; % number of Monte Carlo simulations for each point (value of tuning parameter), used to compute the average value of XT
dBeta0L = 1e4; % product between the difference in propagation constant among consecutive groups and the section length. It can also be regarded as an optimization parameter or fixed it to a "reasonable" value.
tuningParameter_dB_vec = [-80:10:-40 -30:1:-10 -7:3:1 5 10]; % the tuning parameter is the target XT
igsc = 1; % flag to decide whether intra-group strong coupling is enforced via a block-wise random unitary matrix


%% Parameters you likely do not want to change
rndStream = RandStream('mrg32k3a', 'Seed', 0);
nXt = length(tuningParameter_dB_vec);
tuningParamName = "targetXtPerSection_dB";
tuningParamVec = tuningParameter_dB_vec;

nSec = 1; % number of section of the multisectional model. Keep it to 1, since we are looking for the relation between tuning parameter and XT for a single segment
variance_style = "exp_perMode"; % "uniform", "exp_perMode" type of variance
xtMetricsArray = {"SingleValue", "Ferreira"};
if variance_style == "exp_perMode"
    a = -2.25;
    model_parameters = {"a", a}; % parameters of the expm model
else
    model_parameters = [];
end

%% Simulations
xtSingleValueVec = zeros(length(nGroupsVec), nXt);

for nGroupsIdx = 1:length(nGroupsVec)
    if polFlag
        groupSizes = 2*(1:nGroupsVec(nGroupsIdx));
    else
        groupSizes = (1:nGroupsVec(nGroupsIdx));
    end
    nGroups = length(groupSizes);
    nModes = sum(groupSizes);
    startIdx = 1;
    for gg = 1:nGroups
        nModesGroup = groupSizes(gg);
        beta0s(startIdx:startIdx+nModesGroup-1) = dBeta0L*(gg-1);
        startIdx = startIdx + nModesGroup;
    end
    modeIndices = xtHelpers.modeIndicesPGIMMF(nGroups, polFlag); % [1 x nGroups] cell array, where each element is an [1 x nModesPerGroup] array with the indices of the modes which belong to the group
    
    ModelConstructorArgsWihoutTuningParam = {beta0s, dBeta0L, igsc, ...
        polFlag, groupSizes, modeIndices, nSec, variance_style, rndStream};

    [xtSingleValueVec(nGroupsIdx, :), xtMetricsTmp, sectionObj] = expmModel.computeXtVsTuningParameter(@expmModel, tuningParamName, tuningParamVec, ...
        ModelConstructorArgsWihoutTuningParam, xtMetricsArray, nRep, "singleSectionArgs", model_parameters);
    xtMetrics(:, nGroupsIdx, :) = reshape(xtMetricsTmp, length(xtMetricsArray), 1, length(tuningParameter_dB_vec));
end


%% Plots
figure, tiledlayout,
nexttile, plot(tuningParameter_dB_vec, 10*log10(xtSingleValueVec).'), xlabel('target xt inter (dB)'), ylabel('xt inter (dB)'), grid on, lg=legend(string(nGroupsVec)); lg.Title.String = "n groups";
nexttile, plot(tuningParameter_dB_vec, 10*log10(xtSingleValueVec).'-tuningParameter_dB_vec.'), xlabel('target xt inter (dB)'), ylabel('err (dB)'), grid on, lg=legend(string(nGroupsVec)); lg.Title.String = "n groups";

figure, plot(tuningParamVec, 10*log10(xtSingleValueVec)), grid on, hold on,
lgnStr = "Numeric: "+nGroupsVec+" groups";
legend(lgnStr)
xlabel('scaling factor'), ylabel('xt (dB)')

%% Build interpolants
% A set of interpolants is used to numerically compute the relation between
% the tuning parameter, the number of groups and the desired XT

windowSize = 5;
interpolants = cell(1, max(nGroupsVec));
for groupIdx = 1:length(nGroupsVec)
    interpolants{nGroupsVec(groupIdx)} = MultiSectionalModel.buildInterpolant(tuningParamVec, xtSingleValueVec(groupIdx, :), windowSize);
end

%% Test: tuning curve with the generated interpolant
xtVecTest_dB = -80:0.5:0;
tuningParam_test = zeros(1, length(xtVecTest_dB));
testnGroups = 12;
for ii = 1:length(xtVecTest_dB)
    tuningParam_test(ii) = interpolants{testnGroups}(xtVecTest_dB(ii));
end

plot(tuningParam_test, xtVecTest_dB, '--', "DisplayName", "Test generated interpolant"), grid on, hold on
legend

%% SAVE
% Save the numerical relations for later use when simulating a channel.
if igsc == 1
    str = "_with_igsc";
else
    str = "_no_igsc";
end

if dBeta0L == 1
    str2 = "_dBeta0=1";
elseif dBeta0L == 1e4
    str2 = "_dBeta0=1e4";
else
    str2 = "_dBeta0="+round(dBeta0L/1e4)+"e4";
end
    
if polFlag == 1
    str3 = "_with_pol";
else
    str3 = "_no_pol";
end

if variance_style == "uniform"
    str4 = "_unif";
elseif variance_style == "exp_perMode"
    str4 = "_exp_per_mode_a="+a;
end


% save("databases/single_segment_dB_vs_tuning_parameter_expm"+str4+str2+str3+str+".mat", ...
%     "nGroupsVec", "tuningParamName", "tuningParamVec", "xtSingleValueVec", "sectionObj", "xtMetricsArray", "interpolants")

%% Test2: tuning curve with the saved interpolant
if variance_style == "uniform"
    model_parameters_str = [];
elseif variance_style == "exp_perMode"
    model_parameters_str.a = a;
end
expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgsWihoutTuningParam{:}, "model_parameters", model_parameters_str), -20)

xtVecTest_dB = -100:0.5:0;
tuningParam_test = zeros(1, length(xtVecTest_dB));
for ii = 1:length(xtVecTest_dB)
    tuningParam_test(ii) = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgsWihoutTuningParam{:}, "model_parameters", model_parameters_str), xtVecTest_dB(ii));
end

plot(tuningParam_test, xtVecTest_dB, 's--', "DisplayName", "Test loading from the class"), grid on, hold on


%% Test 3: display a power transfer matrix
xtValue_dB = -5; % [dB]
tuningParam = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgsWihoutTuningParam{:}, "model_parameters", model_parameters_str), xtValue_dB);
chObj = expmModel(ModelConstructorArgsWihoutTuningParam{:}, tuningParamName, tuningParam);
xtMetricsNames = {"Ferreira", "SingleValue"};
model_parameters = {}; % parameters of the expm model
concatenationArgs = {"xtMetrics", xtMetricsNames, "nRep", nRep, "singleSectionArgs", model_parameters};
[T_tot, xtMetrics] = chObj.concatenationTransferMatrix(concatenationArgs{:});

figure, imagesc(squeeze(10*log10(mean(abs(T_tot).^2, 1)))), colormap jet; colorbar, axis square