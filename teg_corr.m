function [C, P, DF] = teg_corr(varargin)

% function [C, P, DF] = teg_corr(M[, rank])

M = varargin{1};

if length(varargin) > 1,
    dorank = varargin{2};
else,
    dorank = 0;
end;

m = mean(M, 2);
fn = find(isnan(m));
M(fn, :) = [];

if dorank == 1,
    for iCol = 1:size(M, 2),
        v = M(:, iCol);
        v_order = teg_score_to_order(v);
        M(:, iCol) = v_order;
    end;
end;

try,
    [C, P, DF] = inner_corr(M);
catch,
    C = NaN * ones(size(M, 2), size(M, 2));
    P = ones(size(C));
    DF = zeros(size(C));
end;

function [C, P, DF] = inner_corr(M)
C = NaN * ones(size(M, 2), size(M, 2));
P = NaN * ones(size(M, 2), size(M, 2));
DF = NaN * ones(size(M, 2), size(M, 2));

for iCol1 = 1:size(M, 2),
    for iCol2 = 1:(iCol1 - 1),
        v1 = M(:, iCol1);
        v2 = M(:, iCol2);
        fnn = find(~isnan(v1 + v2));
        v1 = v1(fnn);
        v2 = v2(fnn);
        if isempty(v1),
            continue;
        end;
        v1 = v1 - mean(v1);
        v2 = v2 - mean(v2);
        cov0 = cov(v1, v2);
        sd1 = sqrt(var(v1));
        sd2 = sqrt(var(v2));
        r = cov0(1, 2) / (sd1 * sd2);
        C(iCol1, iCol2) = r;
        N = length(v1);
        DF(iCol1, iCol2) = N - 2;
        if N > 2,
            if abs(r) < 1,
                t = r / sqrt((1 - r^2) / (N - 2));
                p = teg_tpdf(t, N - 2);
            else,
                if isnan(r),
                    p = 1;
                else,
                    p = 0;
                end;
            end;
        else
            p = NaN;
        end;
        P(iCol1, iCol2) = p;
    end;
end;

for iCol1 = 1:size(M, 2),
    for iCol2 = (iCol1 + 1):size(M, 2),
        C(iCol1, iCol2) = C(iCol2, iCol1);
        P(iCol1, iCol2) = P(iCol2, iCol1);
        DF(iCol1, iCol2) = DF(iCol2, iCol1);
    end;
end;

