function [O2L, L, labels] = teg_PCA(varargin)

% function [O2L, L, labels] = teg_PCA(X, varimax, report, labels)

X = varargin{1};
if length(varargin) > 1,
    varimax0 = varargin{2};
else,
    varimax0 = 0;
end;
if length(varargin) > 2,
    report0 = varargin{3};
else,
    report0 = 0;
end;
if length(varargin) > 3,
    itemlabels = varargin{4};
else,
    itemlabels = {};
    for iItem = 1:size(X, 2),
        itemlabels{iItem} = ['obsvar' num2str(iItem)];
    end;
end;

fnn = find(~isnan(mean(X, 2)));
nrows_origX = size(X, 1);
X = X(fnn, :);

X = zscore(X);

[O2L, Lambda] = eig(cov(X));

L = X * O2L;
dL = Lambda;
dL = dL / mean(dL);
fScree = find(dL > 1);

if varimax0 == 1,
    O2L_rot = rotatefactors(O2L(:, fScree));
    L_rot = X * O2L_rot;
    O2L = O2L_rot;
    L = L_rot;
else
    O2L = O2L(:, fScree);
    L = L(:, fScree);
end;

% Flip for interpretation
for iL = 1:size(O2L, 2),
    [ma, ind0] = max(abs(O2L(:, iL)));
    if O2L(ind0, iL) < 0,
        O2L(:, iL) = -O2L(:, iL);
    end;
end;

for iC = 1:size(O2L, 2),
    if report0 == 1,
        fprintf(['Comp ' num2str(iC) ':\n']);
    end;
    for iItem = 1:size(O2L, 1),
        w = O2L(iItem, iC);
        strw = '\t';
        if abs(w) == max(abs(O2L(:, iC))),
            strw = [strw '*'];
            labels{iC} = ['L' num2str(iC) '_' itemlabels{iItem}];
        end;
        if report0 == 1,
            fprintf(['\t' num2str(w, 2) strw '\t' itemlabels{iItem} '\n']);
        end;
    end;
    if report0 == 1,
        fprintf('\n');
    end;
end;

Ltmp = NaN * ones(nrows_origX, size(L, 2));
Ltmp(fnn, :) = L;
