addpath(genpath("./localLibrary"))
addpath(genpath("./databases"))
clearvars

%% User-defined parameters
variance_style = "uniform"; % "uniform", "exp_perMode" type of variance for HelpersN.misc.randomUnitaryMatrix_xt

nGroups = 9; % number of mode groups
polFlag = 1; % 0: no polarizations, 1: with polarizations. It should be kept to 1.
igsc = 1; % flag to add intra-group strong coupling regime via a block-diagonal matrix, whose blocks are uniformly-distributed random unitary matrices
modeIndices = xtHelpers.modeIndicesPGIMMF(nGroups, polFlag);

dBeta0 = 1e4;

xtMetricsNames = {"Ferreira", "SingleValue"}; % XT metrics
nRep = 100; % number of Monte Carlo realizations for the concatenation.
nSec = 100; % number of sections for the multisectional model
targetXtPerSec_dB = -29.5; % [dB] target XT per section

%% Internal parameters
rndStream = RandStream('mrg32k3a', 'Seed', 2);

if variance_style == "exp_perMode"
    a = -2.25;
    model_parameters = {"a", a}; % parameters of the expm model    
    model_parameters_str.a = a;
else
    model_parameters_str = {};
    model_parameters = {};
end

if polFlag
    groupSizes = 2*(1:nGroups); % with polarizations
else
    groupSizes = 1:nGroups; % only spatial modes
end

startIdx = 1;
for gg = 1:nGroups
    nModesGroup = groupSizes(gg);
    beta0s(startIdx:startIdx+nModesGroup-1) = dBeta0*(gg-1);
    startIdx = startIdx + nModesGroup;
end

concatenationArgs = {"xtMetrics", xtMetricsNames, "nRep", nRep, "singleSectionArgs", model_parameters};
ModelConstructorArgs = {beta0s, dBeta0, igsc, ...
        polFlag, groupSizes, modeIndices, nSec, variance_style, rndStream};

tuningParam = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgs{:}, "model_parameters", model_parameters_str), targetXtPerSec_dB);
tuningParamName = "targetXtPerSection_dB";

chObjConcat = expmModel(ModelConstructorArgs{:}, tuningParamName, tuningParam);
[T_tot, xtMetrics] = chObjConcat.concatenationTransferMatrix(concatenationArgs{:});
% T_tot is [nRep, nSec, nModes, nModes]
% xtMetrics is [nRep, nSec, nModes]

%% PLOTS
[~, selIdx] = ismember("SingleValue", string(xtMetricsNames));
xtAvgSingleValueVsDistance_dB = squeeze(10*log10(mean(xtMetrics(selIdx).avgDef1, 1)));

figure, semilogx(xtAvgSingleValueVsDistance_dB), grid on, 
xlabel('Nr of segments'), ylabel('xt avg (dB)'), title('Segment concatenation')

%% Transfer matrix for a specific coupling strength and orientation, and the avg over orientations
xtConcatTarget_dB = -9.52;
nChosenSecs = round(interp1(xtAvgSingleValueVsDistance_dB, 1:length(xtAvgSingleValueVsDistance_dB), xtConcatTarget_dB));
% nChosenSecs = 87;
% nChosenSecs = 1; % TEST

M_single_real = squeeze(T_tot(1,nChosenSecs,:,:));
M_vs_rot = squeeze(T_tot(:,nChosenSecs,:,:));

% Single realization 1
figure, t = tiledlayout;
nexttile, imagesc(abs(M_single_real)), colorbar, title("Single realization: lin scale"), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet;
nexttile, imagesc(10*log10(abs(M_single_real).^2)), colorbar, title("Single realization: log scale"), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet; %clim([-300 0])

% Average
nexttile, imagesc(squeeze(abs(mean(M_vs_rot, 1)))), colorbar, title('Avg over orientation: lin scale'), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet;
nexttile, imagesc(squeeze(10*log10(mean(abs(M_vs_rot).^2, 1)))), colorbar, title('Avg over orientation: log scale'), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet; %clim([-300 0])
title(t, "Coupling matrix")

%% XT: specific coupling strength and orientation
P_mat_dBm = 10*log10(abs(M_single_real).^2);

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};

    if xtMetricName == "Cailabs" || xtMetricName == "Paolo" || xtMetricName == "Ferreira"
        xtInterPerGroupSingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).perGroup(1,nChosenSecs,:,:)));
        xtInterPerGroupAvg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).perGroup(:,nChosenSecs,:,:), 1)));
        plot(xtInterPerGroupSingleReal_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": single real"), grid on, box on, hold on,
        plot(xtInterPerGroupAvg_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": avg over real")
    end

    xtAvg1SingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).avgDef1(1,nChosenSecs,:,:)));
    xtAvg1Avg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef1(:,nChosenSecs,:,:), 1)));

    plot(xtAvg1SingleReal_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 1, single real")
    plot(xtAvg1Avg_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 1, avg over real")

    if xtMetricName == "Cailabs" || xtMetricName == "Paolo"
        xtAvg2SingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).avgDef2(1,nChosenSecs,:,:)));
        xtAvg2Avg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef2(:,nChosenSecs,:,:), 1)));

        plot(xtAvg2SingleReal_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 2, single real")
        plot(xtAvg2Avg_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 2, avg over real")
    end

end
xlabel('Group idx'), ylabel('Inter-group XT (dB)'), legend
title("Single realization. scaling parameter = "+chObjConcat.targetXtPerSection_dB+" dB")

T_tot_end_rep = squeeze(T_tot(end,:,:,:));
T_avg_1 = squeeze(mean(abs(T_tot(:, 1, :, :)).^2, 1));
T_avg_end_seg = squeeze(mean(abs(T_tot(:, nChosenSecs, :, :)).^2, 1));




figure, t=tiledlayout;
nexttile, imagesc(10*log10(abs(squeeze(T_tot_end_rep(1, :,:))).^2)), axis square; colorbar, clim([-80 0]), title('Single real: 1 seg'), colormap jet;
nexttile, imagesc(10*log10(abs(squeeze(T_tot_end_rep(nChosenSecs, :,:))).^2)), axis square; colorbar, clim([-80 0]), title("Single real: "+nChosenSecs+" segs"), colormap jet;
nexttile, imagesc(10*log10(T_avg_1)), axis square; colorbar, clim([-80 0]), title('Mean: 1 seg'), colormap jet;
nexttile, imagesc(10*log10(T_avg_end_seg)), axis square; colorbar, clim([-80 0]), title("Mean: "+nChosenSecs+" segs"), colormap jet;
T_spec_end_seg = abs(squeeze(T_tot_end_rep(nChosenSecs, :, :))).^2;
nexttile, imagesc(10*log10( abs(T_avg_end_seg - T_spec_end_seg) )), axis square; colorbar, clim([-80 0]), title("Difference: single real - mean ("+nChosenSecs+" segs)"), colormap jet;
title(t, "Transfer matrix: "+nChosenSecs+" sections. Each section with scaling factor "+round(tuningParam)+" dB")


%% XT vs distance
colors = lines(10);   % good distinguishable color palette
markers = {'o','s','d','^','v','>','<','p','h','x'};
linestyles = {'-','--',':','-.'};

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    xtAvgVsDistance_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef1, 1)));

    semilogx(1:length(xtAvgVsDistance_dB), xtAvgVsDistance_dB, "DisplayName", xtMetricName+" avg over real", ...
    'Color', colors(xtMetricIdx,:), ...
    'Marker', markers{xtMetricIdx}, ...
    'LineStyle', linestyles{mod(xtMetricIdx-1,length(linestyles))+1}), grid on, hold on
end
xlabel('Section idx'), ylabel('XT (dB)'), legend
title("Concatenation:" + chObjConcat.nSec + " sections")



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Errorbar of average over orientation fixed a coupling strength
%
singleSfIdx = 20;

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    if xtMetricName == "Cailabs" || xtMetricName == "Paolo" || xtMetricName == "Ferreira"

        xtCurr = squeeze(xtMetrics(xtMetricIdx).perGroup(:, singleSfIdx, :));

        % plot(1:nGroups, xtAvgVsDistance_dB, "DisplayName", xtMetricName, ...
        %     'Color', colors(xtMetricIdx,:), ...
        %     'Marker', markers{xtMetricIdx}, ...
        %     'LineStyle', linestyles{mod(xtMetricIdx-1,length(linestyles))+1}), grid on, hold on

        xtHelpers.errorbarLog(1:nGroups, xtCurr.', "dB", "plotOptions", {'x-','DisplayName', xtMetricName}), grid on, box on, hold on,

    end
end
xlabel('Group idx'), ylabel('Inter-group XT (dB)'), legend
title("Errorbars over orientation.")
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% Errorbar of average over orientation for all coupling strengths
%
nn = 10;
figure, t=tiledlayout;
for xtMetricIdx = 1:length(xtMetricsNames)    
    xtMetricName = xtMetricsNames{xtMetricIdx};
    if xtMetricName == "Cailabs" || xtMetricName == "Paolo" || xtMetricName == "Ferreira"
        nexttile
        for secIdx = 1:nn:nSec
            xtHelpers.errorbarLog(1:nGroups, squeeze(xtMetrics(xtMetricIdx).perGroup(:, secIdx, :)).', "dB", "plotOptions", {'x-','DisplayName', xtMetricName}), grid on, box on, hold on,
        end
        xlabel('Group idx'), ylabel('Inter-group XT (dB)'), title("Cailabs def")
    end
end
xlabel('Group idx'), ylabel('Inter-group XT (dB)'), legend
title(t, "Errorbars over orientation, for different coupling strengths")
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %












%% Single section, sweep over xt
nSec = 1;
targetXtPerSecVec_dB = -50:5:0;
nSfVec = length(targetXtPerSecVec_dB);

xtMetricsArray = cell(1, nSfVec); % cell of cells
for sfIdx = 1:nSfVec
    targetXtPerSec_dB = targetXtPerSecVec_dB(sfIdx); % [dB]
    tuningParam = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgs{:}, "model_parameters", model_parameters_str), targetXtPerSec_dB);
    tuningParamName = "targetXtPerSection_dB";    

    ModelConstructorArgs = {beta0s, dBeta0, igsc, ...
        polFlag, groupSizes, modeIndices, nSec, variance_style, rndStream};


    chObj1SecSweep = expmModel(ModelConstructorArgs{:}, tuningParamName, tuningParam);
    sweepArgs = {"xtMetrics", xtMetricsNames, "nRep", nRep, "singleSectionArgs", model_parameters};
    [T_tot, xtMetricsArray{sfIdx}] = chObj1SecSweep.concatenationTransferMatrix(sweepArgs{:});
    % T_tot is [nRep, nSec, nModes, nModes]
    % xtMetricsArray is [nSfVec, nXtDefs] of structures
end


%% XT vs coupling strength
colors = lines(10);   % good distinguishable color palette
markers = {'o','s','d','^','v','>','<','p','h','x'};
linestyles = {'-','--',':','-.'};


figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    
    xtAvgVsDistance_dB = zeros(1, length(xtMetricsArray));
    for sfIdx = 1:length(xtMetricsArray)
        xtAvgVsDistance_dB(sfIdx) = squeeze(10*log10(mean(xtMetricsArray{sfIdx}(xtMetricIdx).avgDef1, 1)));
    end

    plot(targetXtPerSecVec_dB, xtAvgVsDistance_dB, "DisplayName", xtMetricName+" avg over real", ...
    'Color', colors(xtMetricIdx,:), ...
    'Marker', markers{xtMetricIdx}, ...
    'LineStyle', linestyles{mod(xtMetricIdx-1,length(linestyles))+1}), grid on, hold on
end
plot(targetXtPerSecVec_dB, targetXtPerSecVec_dB, "--", "DisplayName", "Ideal")
xlabel('Target XT'), ylabel('XT (dB)'), legend



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Errorbar of average over orientation fixed a coupling strength
%
singleSfIdx = 4;

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    if xtMetricName == "Cailabs" || xtMetricName == "Paolo" || xtMetricName == "Ferreira"

        xtCurr = squeeze(xtMetricsArray{singleSfIdx}(xtMetricIdx).perGroup);

        % plot(1:nGroups, xtAvgVsDistance_dB, "DisplayName", xtMetricName, ...
        %     'Color', colors(xtMetricIdx,:), ...
        %     'Marker', markers{xtMetricIdx}, ...
        %     'LineStyle', linestyles{mod(xtMetricIdx-1,length(linestyles))+1}), grid on, hold on

        xtHelpers.errorbarLog(1:nGroups, xtCurr.', "dB", "plotOptions", {'x-','DisplayName', xtMetricName}), grid on, box on, hold on,

    end
end
xlabel('Group idx'), ylabel('Inter-group XT (dB)'), legend
title("Errorbars over orientation, target XT = "+targetXtPerSecVec_dB(singleSfIdx)+" dB")
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% Errorbar of average over orientation for all coupling strengths
%
nn = 1;
figure, t=tiledlayout;
for xtMetricIdx = 1:length(xtMetricsNames)    
    xtMetricName = xtMetricsNames{xtMetricIdx};
    if xtMetricName == "Cailabs" || xtMetricName == "Paolo" || xtMetricName == "Ferreira"
        nexttile
        for tuningParamIdx = 1:nn:nSfVec
            xtHelpers.errorbarLog(1:nGroups, squeeze(xtMetricsArray{tuningParamIdx}(xtMetricIdx).perGroup).', "dB", "plotOptions", {'x-','DisplayName', xtMetricName}), grid on, box on, hold on,
        end
        xlabel('Group idx'), ylabel('Inter-group XT (dB)'), title("Cailabs def")
    end
end
xlabel('Group idx'), ylabel('Inter-group XT (dB)'), legend
title(t, "Errorbars over orientation, for different coupling strengths")
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


