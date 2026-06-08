classdef xtHelpers
%XTHELPERS  Static helper utilities for mode-index bookkeeping, XT plotting, and random matrices.
%
%   xtHelpers is a collection of static methods used throughout the
%   multi-sectional fibre model framework. All methods are static;
%   instantiation is not required.
%
%   Methods:
%     modeIndicesPGIMMF    - generate per-group mode indices for a parabolic GI-MMF
%     drawGroupLimits      - overlay mode-group boundaries on an imagesc() plot
%     errorbarLog          - errorbar plot on a log10 or dB scale
%     randomUnitaryMatrix  - Haar-distributed random unitary matrix
%
%   Example:
%     modeIdx = xtHelpers.modeIndicesPGIMMF(4, true);
%     imagesc(abs(H).^2);
%     xtHelpers.drawGroupLimits(modeIdx);
%
%   Author:      Paolo Carniello
%   Affiliation: Technical University of Munich
%   Date:        2026-06-08
%
%   See also MULTISECTIONALMODEL, EXPMMODEL.

    methods (Static)
        function modeIndices = modeIndicesPGIMMF(nGroups, polFlag)
            %MODEINDICES_PGIMMF  Per-group mode indices of a parabolic graded-index MMF.
            %
            %   MODEINDICES = MODEINDICES_PGIMMF(NGROUPS, POLFLAG) returns the mode
            %   index cell array for a parabolic graded-index multimode fibre (PGI-MMF).
            %   Group g contains g modes (or 2*g modes when polarisation is included).
            %   Indices are assigned sequentially starting from 1.
            %
            %   Inputs:
            %     NGROUPS - [1 x 1 double] number of mode groups
            %     POLFLAG - [1 x 1 logical] if true, group g has 2*g modes;
            %               otherwise g modes
            %
            %   Output:
            %     MODEINDICES - [1 x nGroups cell] each element is a
            %                   [1 x nModesPerGroup double] vector of global mode
            %                   indices belonging to that group
            %
            %   Example:
            %     modeIdx = xtHelpers.modeIndicesPGIMMF(3, false);
            %     % modeIdx = {[1], [2 3], [4 5 6]}
            %
            %   See also EXPMMODEL, MULTISECTIONALMODEL.
            cnt = 1;
            if polFlag
                nModesPerGroup = 2*(1:nGroups);
            else
                nModesPerGroup = 1:nGroups;
            end
            modeIndices = cell(1, nGroups);
            for groupIdx = 1:nGroups
                modeIndices{groupIdx} = cnt:(cnt+nModesPerGroup(groupIdx)-1);
                cnt = cnt + nModesPerGroup(groupIdx);
            end
        end


        function drawGroupLimits(modeIndices, varargin)
            %DRAWGROUPLIMITS  Overlay mode-group boundary lines on an imagesc() plot.
            %
            %   DRAWGROUPLIMITS(MODEINDICES) draws white horizontal and vertical lines
            %   on the current axes at the boundaries between mode groups, suitable
            %   for annotating an imagesc() plot of an [nModes x nModes] matrix.
            %
            %   DRAWGROUPLIMITS(MODEINDICES, 'LineWidth', W) uses line width W.
            %
            %   Input:
            %     MODEINDICES - [1 x nGroups cell] each element is a
            %                   [1 x nModesPerGroup double] mode index vector
            %
            %   Name-Value Arguments:
            %     'LineWidth' - [1 x 1 double] width of the boundary lines;
            %                   default 3
            %
            %   Example:
            %     imagesc(10*log10(abs(H).^2));
            %     xtHelpers.drawGroupLimits(modeIndices, 'LineWidth', 2);
            %
            %   See also MODEINDICES_PGIMMF, IMAGESC.
            p = inputParser();
            p.addParameter("LineWidth", 3);
            p.parse(varargin{:});



            nGroups = length(modeIndices);
            idx = zeros(1, nGroups);
            for g = 1:nGroups
                idx(g) = modeIndices{g}(end);
            end

            nModes = max(modeIndices{end});

            for k = 1:length(idx)
                x = idx(k) + 0.5;

                % Vertical line (between columns)
                line([x x], [0.5 nModes+0.5], 'Color', 'w', 'LineWidth', p.Results.LineWidth)

                % Horizontal line (between rows)
                line([0.5 nModes+0.5], [x x], 'Color', 'w', 'LineWidth', p.Results.LineWidth)
            end
        end

        function [meanVal_log, dX_maxmin_log, dX_stddev_log] = errorbarLog(xAxis, xtPerGroup, unit, varargin)
            %ERRORBARLOG  Errorbar plot of crosstalk on a log10 or dB scale.
            %
            %   [MEANVAL_LOG, DX_MAXMIN_LOG, DX_STDDEV_LOG] = ...
            %       ERRORBARLOG(XAXIS, XTPERGROUP, UNIT) plots an errorbar of
            %       XTPERGROUP versus XAXIS on a logarithmic scale and returns the
            %       plotted statistics.
            %
            %   Error bars are computed in the log domain and the bar length is either 
            %   equal to log10(max)-log10(min), where max and min are the maximum and 
            %   minimum values of the XT realizations, or to the stddev of the XT realizations.
            %
            %   [...] = ERRORBARLOG(..., Name, Value) accepts optional arguments:
            %
            %   Name-Value Arguments:
            %     'barLength'  - [string] "maxmin" uses the half max-min range;
            %                    "stddev" uses mean+std - mean in log scale;
            %                    default "maxmin"
            %     'plotOptions'- [cell] extra arguments forwarded to errorbar()
            %
            %   Inputs:
            %     XAXIS      - [1 x nxPoints double] x-axis values
            %     XTPERGROUP - [nxPoints x nRep double] crosstalk samples in linear
            %                  scale (not dB); each row is one x-axis point,
            %                  each column one Monte-Carlo realisation
            %     UNIT       - [string] "dB" multiplies log10 outputs by 10;
            %                  any other value returns plain log10
            %
            %   Outputs:
            %     MEANVAL_LOG   - [nxPoints x 1 double] mean XT in log10 or dB
            %     DX_MAXMIN_LOG - [nxPoints x 1 double] half max-min bar length
            %                     in log10 or dB
            %     DX_STDDEV_LOG - [nxPoints x 1 double] std-dev bar half-length
            %                     in log10 or dB
            %
            %   Example:
            %     xtHelpers.errorbarLog(nGroupsVec, xtData, "dB", ...
            %         'barLength', 'stddev', 'plotOptions', {'b-o'});
            %
            %   See also ERRORBAR, DRAWGROUPLIMITS.
            p = inputParser();
            p.addParameter("barLength", "maxmin") % "maxmin" or "stddev". "maxmin" sets as errorbar length the difference between max and min values (in log10 or dB), stddev the difference between mean+std and mean (in log10 or dB)
            p.addParameter("plotOptions", []); % things to be passed to errorbar as {:}
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            meanVal = mean(xtPerGroup, 2);
            meanVal_log = log10(meanVal);
            dX_maxmin_log = (log10(max(xtPerGroup, [], 2))-log10(min(xtPerGroup, [], 2)))/2;

            stddev = sqrt(var(xtPerGroup, [], 2));
            dX_stddev_log = log10(meanVal+stddev) - log10(meanVal);

            if unit == "dB"
                meanVal_log = 10*meanVal_log;
                dX_maxmin_log = 10*dX_maxmin_log;
                dX_stddev_log = 10*dX_stddev_log;
            end
            if p.Results.barLength == "maxmin"
                errorbar(xAxis, meanVal_log, dX_maxmin_log, p.Results.plotOptions{:})
            elseif p.Results.barLength == "stddev"
                errorbar(xAxis, meanVal_log, dX_stddev_log, p.Results.plotOptions{:})
            end
        end

        function unitaryMat = randomUnitaryMatrix(N, rndStream)
            %RANDOMUNIATRYMATRIX  Generate an N x N Haar-random unitary matrix.
            %
            %   UNITARYMAT = RANDOMUNIATRYMATRIX(N, RNDSTREAM) returns an [N x N]
            %   unitary matrix drawn uniformly from U(N) according to the Haar measure,
            %   using the numerically stable QR-based algorithm from [Mezzadri, 2007].
            %
            %   Inputs:
            %     N         - [1 x 1 double] dimension of the desired square matrix
            %     RNDSTREAM - [RandStream] random stream for the Gaussian draws
            %
            %   Output:
            %     UNITARYMAT - [N x N double] Haar-random unitary matrix
            %
            %   Reference:
            %     F. Mezzadri, "How to Generate Random Matrices from the Classical
            %     Compact Groups", Notices of the AMS, 54(5), pp. 592-604, 2007.
            %
            %   Example:
            %     rng = RandStream('mrg32k3a', 'Seed', 42);
            %     U   = xtHelpers.randomUnitaryMatrix(6, rng);
            %     norm(U * U' - eye(6))   % should be ~eps
            %
            %   See also INTRAGROUPSTRONGCOUPLING.

            gaussMat = (randn(rndStream, N, N) + 1i*randn(rndStream, N, N))/sqrt(2);

            [unitaryMat, r] = qr(gaussMat);
            lambda = diag(r); lambda = diag(lambda./abs(lambda));
            unitaryMat = unitaryMat * lambda;
        end
    end
end
