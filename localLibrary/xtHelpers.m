classdef xtHelpers
    methods (Static)
        function modeIndices = modeIndicesPGIMMF(nGroups, polFlag)
            % Mode indices of a parabolic graded-index multimode fiber, in
            % the format used in this class.

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
            % Draw white boundaries on an imagesc() relative to an nModes x nModes
            % matrix relative to a fiber with group sizes from 1 to
            % nGroups.

            % INPUTS:
            % -modeIndices: [1 x nGroups] cell array, with [1 x nModesPerGroup] elements

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
            % INPUT:
            % xAxis: [1 x nxPoints]
            % xtPerGroup: [nxPoints x nRep]
            % The bar stretches to log10(mean(xtPerGroup))+dX, where
            % log10(mean(xtPerGroup))+dX is such that
            % log10(mean(xtPerGroup))+dX = log10(mean(xtPerGroup)+d), where d =
            % (max(xtPerGroup)-min(xtPerGroup))/2. Hence, dX =
            % log10(mean(xtPerGroup)+d) - log10(mean(xtPerGroup));

            p = inputParser();
            p.addParameter("barLength", "maxmin") % "maxmin" or "stddev". "maxmin" sets as errorbar length the difference between max and min values (in log10 or dB), stddev the difference between mean+std and mean (in log10 or dB)
            p.addParameter("plotOptions", []); % things to be passed to errorbar as {:}
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            meanVal = mean(xtPerGroup, 2);
            meanVal_log = log10(meanVal);
            % d = (max(xtPerGroup, [], 2)-min(xtPerGroup, [], 2))/2;
            % dX_log = log10(mean(xtPerGroup, 2)+d) - meanVal_log;
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
            % INPUTS:
            % N: dimension of the desired [N x N] square matrix
            % variance: see the code. Not sure it actually makes a difference. My
            % desire is to have unitary matrix which couple more and less based on a
            % certain parameter, but I'm not sure which is such parameter.

            % Algorithm from pag.597 of Mezzadri2007, ensuring stability
            
            gaussMat = (randn(rndStream, N, N) + 1i*randn(rndStream, N, N))/sqrt(2);

            [unitaryMat, r] = qr(gaussMat);
            lambda = diag(r); lambda = diag(lambda./abs(lambda));
            unitaryMat = unitaryMat * lambda;
        end
    end
end