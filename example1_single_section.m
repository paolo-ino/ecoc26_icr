% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Example script showing how to use the library to generate the transfer
% matrix for linear mode coupling of a concatenation of section of a graded
% index multimode fiber.
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

%% Parameters you might want to change
targetXtPerSec_dB = -20; % [dB] desired value of average fiber crosstalk as in [Carniello26, p.2]
nGroups = 9; % number of mode groups for a graded-index multimode fiber
polFlag = 1; % flag to decide if polarizations are considered as well, or only spatial modes
nRep = 100; % number of Monte Carlo simulations for each point (value of tuning parameter), used to compute the average value of XT
dBeta0L = 1e4; % product between the difference in propagation constant among consecutive groups and the section length. It can also be regarded as an optimization parameter or fixed it to a "reasonable" value.
igsc = 1; % flag to decide whether intra-group strong coupling is enforced via a block-wise random unitary matrix
variance_style = "uniform"; % "uniform", "exp_perMode" type of variance

%% Parameters you likely do not want to change
rndStream = RandStream('mrg32k3a', 'Seed', 0);
nSec = 1; % number of section of the multisectional model. Keep it to 1, since we are looking for the relation between tuning parameter and XT for a single segment
xtMetricsArray = {"SingleValue", "Ferreira"};
if variance_style == "exp_perMode"
    a = -2.25;
    model_parameters = {"a", a}; % parameters of the expm model
else
    model_parameters = [];
end

if polFlag
    groupSizes = 2*(1:nGroups);
else
    groupSizes = (1:nGroups);
end
startIdx = 1;
for gg = 1:nGroups
    nModesGroup = groupSizes(gg);
    beta0sL(startIdx:startIdx+nModesGroup-1) = dBeta0L*(gg-1);
    startIdx = startIdx + nModesGroup;
end
modeIndices = xtHelpers.modeIndicesPGIMMF(nGroups, polFlag); % [1 x nGroups] cell array, where each element is an [1 x nModesPerGroup] array with the indices of the modes which belong to the group
xtMetricsNames = {"Ferreira", "SingleValue"};

concatenationArgs = {"xtMetrics", xtMetricsNames, "nRep", nRep, "singleSectionArgs", model_parameters};
ModelConstructorArgs = {beta0sL, dBeta0L, igsc, ...
        polFlag, groupSizes, modeIndices, nSec, variance_style, rndStream};

% Obtain the tuning parameter of the model from the target XT
tuningParam = expmModel.getTuningParameterFromXt(expmModel(ModelConstructorArgs{:}, "model_parameters", model_parameters_str), targetXtPerSec_dB);
tuningParamName = "targetXtPerSection_dB";

%% Simulation
chObjConcat = expmModel(ModelConstructorArgs{:}, tuningParamName, tuningParam);
[T_tot, xtMetrics] = chObjConcat.concatenationTransferMatrix(concatenationArgs{:});

%% Plots
chosenSecIdx = nSec; % section index for which plotting the various quantities
P_mat_mean_dB = squeeze(10*log10(mean(abs(T_tot(:, chosenSecIdx, :, :)).^2, 1)));

% Transfer matrix
figure, imagesc(P_mat_mean_dB), colormap jet; colorbar, axis square
a=colorbar; ylabel(a,'Relative power (dB)','FontSize',16,'Rotation',270);
xtHelpers.drawGroupLimits(modeIndices);
xlabel('Input mode index'), ylabel('Output mode index')
title("Average power transfer matrix")

% Crosstalk per group
chosenRealIdx = 1; % realization index for which plotting the various quantities

figure,
for xtMetricIdx = 1:length(xtMetricsNames)
    xtMetricName = xtMetricsNames{xtMetricIdx};

    if xtMetricName == "Ferreira"
        xtInterPerGroupSingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).perGroup(chosenRealIdx,chosenSecIdx,:,:)));
        xtInterPerGroupAvg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).perGroup(:,chosenSecIdx,:,:), 1)));
        plot(xtInterPerGroupSingleReal_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": single real"), grid on, box on, hold on,
        plot(xtInterPerGroupAvg_dB, 'x-','DisplayName', "Per group "+ xtMetricName+": avg over real")
    end

    xtAvg1SingleReal_dB = squeeze(10*log10(xtMetrics(xtMetricIdx).avgDef1(chosenRealIdx,chosenSecIdx,:,:)));
    xtAvg1Avg_dB = squeeze(10*log10(mean(xtMetrics(xtMetricIdx).avgDef1(:,chosenSecIdx,:,:), 1)));

    plot(xtAvg1SingleReal_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 1, single real")
    plot(xtAvg1Avg_dB*ones(1, nGroups), '--','DisplayName', xtMetricName + " avg per fiber 1, avg over real")
end
xlabel('Group idx'), ylabel('XT (dB)'), legend