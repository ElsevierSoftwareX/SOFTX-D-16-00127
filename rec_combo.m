function [meanvec, idmat, Nvec, SDvec, levels, levelVals] = rec_combo(varargin)

% function [meanvec, idmat, Nvec, SDvec, levels, levelVals] = rec_combo(depvec, indepmat, useMedian, removeOutliers, removeTrendPower)

depvec = varargin{1};
indepvec = varargin{2};
if length(varargin) > 2,
    useMedian = varargin{3}; % Or alt
else,
    useMedian = 0;
end;
if length(varargin) > 3,
    removeOutliers = varargin{4};
else,
    removeOutliers = 0;
end;
if length(varargin) > 4,
    removeTrend = varargin{5};
else,
    removeTrend = 0;
end;

% Make sure all combinations are included and estimated if necessary
levels = [];
levelVals = [];
for iVar = 1:size(indepvec, 2),
    u = unique(indepvec(:, iVar));
    levels = [levels length(u)];
    for iu = 1:length(u),
        levelVals = [levelVals u(iu)];
    end;
end;

indepvec0 = rec_comb(levels, 0);
for iVar = 1:size(indepvec0, 2),
    u = unique(indepvec(:, iVar));
    tmp = indepvec0(:, iVar);
    for iu = 1:length(u),
        f = find(tmp == iu);
        indepvec0(f, iVar) = u(iu);
    end;
end;
depvec0 = NaN * ones(size(indepvec0, 1), 1);
depvec = [depvec0; depvec];
indepvec = [indepvec0; indepvec];
[meanvec, idmat, Nvec, SDvec] = rec_combo_inner(depvec, indepvec, [], 1, [], [], [], [], useMedian, removeOutliers, removeTrend);

function [meanvec, idmat, Nvec, SDvec] = rec_combo_inner(depvec, indepvec, meanvec, depth, idmat, idval, Nvec, SDvec, useMedian, removeOutliers, removeTrend)

if (depth == size(indepvec, 2)),
    u = unique(indepvec(:, depth));
    for iu = 1:length(u),
        if isnan(u(iu)),
            continue;
        end;
        f = find(indepvec(:, depth) == u(iu));
        idval(1, depth) = iu;
        selected = depvec(f);
        selected(find(isnan(selected))) = [];
        if removeOutliers == 1,
            foutlier = inner_remove_outliers(selected);
            selected(foutlier) = [];
        end;
        if removeTrend > 0,
            % For using SD
            msel = mean(selected);
            ToT = ones(length(selected), 1);
            t = 1:length(ToT);
            t = t(:);
            t = t ./ max(t);
            t = t - mean(t);
            for p = 1:removeTrend,
                ToT = [ToT t .^ p];
            end;
            try,
                b = inv(ToT'*ToT)*ToT'*selected;
                if length(find(isnan(b))) > 0,
                    fprintf('No trend removal: NaNs\n');
                    dfdfd=0;
                end;
                pred = ToT * b;
                selected = msel + selected - pred;                
            catch,
                fprintf('No trend removal: error\n');
                dfdf=0;
            end;
        end;
        if ~isempty(selected),
            if useMedian == 1,
                newmeanvec = median(selected(~isnan(selected)));
            elseif useMedian == 2, % ex-Gaussian
                params = egfit(selected);
                newmeanvec = params(1);
            elseif useMedian == 3, % speed
                newmeanvec = mean(1 ./ selected(~isnan(selected)));
            else,
                newmeanvec = mean(selected(~isnan(selected)));
            end;
            meanvec = [meanvec newmeanvec];
            SDvec = [SDvec sqrt(var(selected(~isnan(selected))))];
        else,
            meanvec = [meanvec NaN];
            SDvec = [SDvec NaN];
        end;
        idmat = [idmat; idval];
        Nvec = [Nvec length(selected)];
    end;
else,
    u = unique(indepvec(:, depth));
    for iu = 1:length(u),
        if isnan(u(iu)),
            continue;
        end;
        f = find(indepvec(:, depth) == u(iu));
        idval(1, depth) = iu;
        [meanvec, idmat, Nvec, SDvec] = rec_combo_inner(depvec(f), indepvec(f, :), meanvec, depth + 1, idmat, idval, Nvec, SDvec, useMedian, removeOutliers, removeTrend);
    end;
end;

function foutlier = inner_remove_outliers(selected)
z = zscore(selected);
foutlier = find(abs(z) > 3);
