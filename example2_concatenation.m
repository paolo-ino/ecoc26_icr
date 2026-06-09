% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Example script showing how to use the library to generate the transfer
% matrix for linear mode coupling of a concatenation of section of a graded
% index multimode fiber. See readme.md and [Carniello26] for more details
% on the models.
% 
% The current choice of parameters allows to reproduce the results for the
% proposed expm models of [Carniello26, Fig.2-3].
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


addpath(genpath("./localLibrary"))
addpath(genpath("./databases"))
clearvars

%% Parameters you might want to change
variance_style = "uniform"; % "uniform", "exp_perMode" type of variance for HelpersN.misc.randomUnitaryMatrix_xt

nGroups = 9; % number of mode groups
polFlag = 1; % 0: no polarizations, 1: with polarizations. It should be kept to 1.
igsc = 1; % flag to add intra-group strong coupling regime via a block-diagonal matrix, whose blocks are uniformly-distributed random unitary matrices
modeIndices = xtHelpers.modeIndicesPGIMMF(nGroups, polFlag);

dBeta0L = 1e4; % product between the difference in propagation constant among consecutive groups and the section length. It can also be regarded as an optimization parameter or fixed it to a "reasonable" value.

xtMetricsNames = {"Ferreira", "SingleValue"}; % XT metrics
nRep = 100; % number of Monte Carlo realizations for the concatenation.
nSec = 100; % number of sections for the multisectional model
targetXtPerSec_dB = -29.5; % [dB] target XT per section

%% Parameters you probably don't want to change
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
    beta0s(startIdx:startIdx+nModesGroup-1) = dBeta0L*(gg-1);
    startIdx = startIdx + nModesGroup;
end

concatenationArgs = {"xtMetrics", xtMetricsNames, "nRep", nRep, "singleSectionArgs", model_parameters};
ModelConstructorArgs = {beta0s, dBeta0L, igsc, ...
        polFlag, groupSizes, modeIndices, nSec, variance_style, rndStream};

tuningParam = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgs{:}, "model_parameters", model_parameters_str), targetXtPerSec_dB);
tuningParamName = "targetXtPerSection_dB";

%% Simulations
chObjConcat = expmModel(ModelConstructorArgs{:}, tuningParamName, tuningParam);
[T_tot, xtMetrics] = chObjConcat.concatenationTransferMatrix(concatenationArgs{:});
% T_tot is [nRep, nSec, nModes, nModes]
% xtMetrics is [nRep, nSec, nModes]

%% Plots: Transfer matrices
[~, selIdx] = ismember("SingleValue", string(xtMetricsNames));
xtAvgSingleValueVsDistance_dB = squeeze(10*log10(mean(xtMetrics(selIdx).avgDef1, 1)));

figure, semilogx(xtAvgSingleValueVsDistance_dB), grid on, 
xlabel('Nr of segments'), ylabel('xt avg (dB)'), title('Segment concatenation')

%% Plots: Transfer matrix for a specific coupling strength and realization, and the average over realizations
xtConcatTarget_dB = -9.52;
nChosenSecs = round(interp1(xtAvgSingleValueVsDistance_dB, 1:length(xtAvgSingleValueVsDistance_dB), xtConcatTarget_dB));

M_single_real = squeeze(T_tot(1,nChosenSecs,:,:));
M_vs_rot = squeeze(T_tot(:,nChosenSecs,:,:));

% Single realization
figure, t = tiledlayout;
nexttile, imagesc(abs(M_single_real)), colorbar, title("Single realization: lin scale"), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet;
nexttile, imagesc(10*log10(abs(M_single_real).^2)), colorbar, title("Single realization: log scale"), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet; %clim([-300 0])

% Average
nexttile, imagesc(squeeze(abs(mean(M_vs_rot, 1)))), colorbar, title('Avg over realizations: lin scale'), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet;
nexttile, imagesc(squeeze(10*log10(mean(abs(M_vs_rot).^2, 1)))), colorbar, title('Avg over realizations: log scale'), xtHelpers.drawGroupLimits(modeIndices); axis square; colormap jet; %clim([-300 0])
title(t, "Coupling matrix")

%% Plots: Crosstalk per group
P_mat_dBm = 10*log10(abs(M_single_real).^2);

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};

    if xtMetricName == "Ferreira"
        xtInterPerGroupSingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).perGroup(1,nChosenSecs,:)));
        xtInterPerGroupAvg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).perGroup(:,nChosenSecs,:), 1)));
        plot(xtInterPerGroupSingleReal_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": single real"), grid on, box on, hold on,
        plot(xtInterPerGroupAvg_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": avg over realizations")
    end

    xtAvg1SingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).avgDef1(1,nChosenSecs)));
    xtAvg1Avg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef1(:,nChosenSecs), 1)));

    plot(xtAvg1SingleReal_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber, single real")
    plot(xtAvg1Avg_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber, avg over realizations")

end
xlabel('Group idx'), ylabel('XT (dB)'), legend
title("Single realization. XT per section = "+ round(chObjConcat.targetXtPerSection_dB*100)/100+" dB")

T_tot_end_rep = squeeze(T_tot(end,:,:,:));
T_avg_1 = squeeze(mean(abs(T_tot(:, 1, :, :)).^2, 1));
T_avg_end_seg = squeeze(mean(abs(T_tot(:, nChosenSecs, :, :)).^2, 1));




figure, t=tiledlayout;
nexttile, imagesc(10*log10(abs(squeeze(T_tot_end_rep(1, :,:))).^2)), axis square; colorbar, clim([-80 0]), title('Single real: 1 seg'), colormap jet;
nexttile, imagesc(10*log10(abs(squeeze(T_tot_end_rep(nChosenSecs, :,:))).^2)), axis square; colorbar, clim([-80 0]), title("Single real: "+nChosenSecs+" segs"), colormap jet;
nexttile, imagesc(10*log10(T_avg_1)), axis square; colorbar, clim([-80 0]), title('Mean: 1 seg'), colormap jet;
nexttile, imagesc(10*log10(T_avg_end_seg)), axis square; colorbar, clim([-80 0]), title("Mean: "+nChosenSecs+" segs"), colormap jet;
T_spec_end_seg = abs(squeeze(T_tot_end_rep(nChosenSecs, :, :))).^2;
title(t, "Transfer matrix: "+nChosenSecs+" sections. Each section with scaling factor "+round(tuningParam)+" dB")


%% Plots: Crosstalk vs distance
colors = lines(10);   % good distinguishable color palette
markers = {'o','s','d','^','v','>','<','p','h','x'};
linestyles = {'-','--',':','-.'};

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    xtAvgVsDistance_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef1, 1)));

    semilogx(1:length(xtAvgVsDistance_dB), xtAvgVsDistance_dB, "DisplayName", xtMetricName+" avg over realizations", ...
    'Color', colors(xtMetricIdx,:), ...
    'Marker', markers{xtMetricIdx}, ...
    'LineStyle', linestyles{mod(xtMetricIdx-1,length(linestyles))+1}), grid on, hold on
end
xlabel('Section idx'), ylabel('XT (dB)'), legend
title("Concatenation: " + chObjConcat.nSec + " sections")


%% Plots: Errorbar of average over realizations
%
figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};
    if xtMetricName == "Ferreira"

        xtCurr = squeeze(xtMetrics(xtMetricIdx).perGroup(:, nChosenSecs, :));
        xtHelpers.errorbarLog(1:nGroups, xtCurr.', "dB", "plotOptions", {"x-","DisplayName", sprintf("Per group " +xtMetricName)}), grid on, box on, hold on,

    end

    if xtMetricName == "SingleValue"
        xtHelpers.errorbarLog(1:nGroups, repmat(xtMetrics(xtMetricIdx).avgDef1(:, nChosenSecs).', nGroups, 1), "dB", "plotOptions", {"x-","DisplayName", sprintf("Avg per fiber"+xtMetricName)}), grid on, box on, hold on,
    end
end
xlabel('Group idx'), ylabel('XT (dB)'), legend
title("Errorbars over realizations")
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

